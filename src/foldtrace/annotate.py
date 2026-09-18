"""The fetch layer: build the annotation table the detectors consume.

Everything else in foldtrace is offline by design, which is what lets the same
analysis run on a laptop, in CI, and in a reviewer's hands from an archived
deposit. That property is worth keeping, so the network lives here, at the edge,
and nowhere else.

Three commitments, each of which exists because of a way this step goes wrong:

*Provenance, not silence.* Every requested accession gets a row. A record that
could not be fetched is marked with why -- ``not_found``, ``http_error``,
``withdrawn`` -- and never silently omitted, because an absent row reads
downstream as an absent feature, and a failed lookup is not a biological
absence. The count of rows out always equals the count of accessions in.

*Withdrawn records are recovered, not dropped.* Entries removed from active
UniProtKB still carry organism and protein-name information in UniParc, and
they are exactly the under-annotated proteins a structure-first search is for:
one published organophosphate-hydrolase set had 1,435 of them.

*A manifest the freeze can hash.* The register locks the sites file and the
prior-art claims; it should also lock where the annotations came from. This
writes a manifest recording each endpoint, the query date, the counts by
status, and the SHA256 of the output, so ``prime freeze --also manifest.json``
puts the inputs under the same seal as everything else.

Network use is confined to :func:`fetch_json`. Every parser is pure and tested
against recorded payloads. Run ``foldtrace annotate --probe ACCESSION`` first:
it fetches one record and prints the raw JSON, so the response shape can be
confirmed against a real answer before a thousand requests are made on the
strength of an assumed one.

Stdlib only.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone

UNIPROT = "https://rest.uniprot.org/uniprotkb/{acc}.json"
UNIPARC = "https://rest.uniprot.org/uniparc/search?query=upi:{acc}&format=json"
UNIPARC_BY_ACC = "https://rest.uniprot.org/uniparc/search?query={acc}&format=json"
INTERPRO_PFAM = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam}/"

ANNOTATION_COLUMNS = ["target", "name", "pfam", "clan", "organism", "taxon_id",
                      "length", "status", "source"]

USER_AGENT = "foldtrace-annotate/0.1 (+https://github.com/headsparkle/foldtrace)"


class AnnotateError(Exception):
    """Raised for input problems; network failures become row statuses instead."""


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# --------------------------------------------------------------------- network
def fetch_json(url: str, timeout: float = 30.0, retries: int = 2,
               pause: float = 0.34) -> tuple[dict | list | None, str]:
    """The only function here that touches the network.

    Returns ``(payload, status)``. Status is ``ok``, ``not_found`` for a 404,
    or ``http_error:<code>`` / ``network_error`` otherwise -- never an
    exception, because a failed lookup is a recorded outcome rather than a
    crash. ``pause`` keeps the default request rate polite; UniProt and
    InterPro both ask for it.
    """
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                               "Accept": "application/json"})
    last = "network_error"
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
            time.sleep(pause)
            return payload, "ok"
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None, "not_found"
            last = f"http_error:{exc.code}"
            if exc.code < 500:
                break           # a client error will not fix itself
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = f"network_error:{type(exc).__name__}"
        if attempt < retries:
            time.sleep(pause * (2 ** attempt))
    return None, last


# ---------------------------------------------------------------------- cache
class Cache:
    """One JSON file per URL. Makes a re-run free and an audit possible."""

    def __init__(self, directory: str | None):
        self.dir = directory
        if directory:
            os.makedirs(directory, exist_ok=True)

    def _path(self, url: str) -> str:
        return os.path.join(self.dir, hashlib.sha256(url.encode()).hexdigest()[:32] + ".json")

    def get(self, url: str):
        if not self.dir:
            return None
        p = self._path(url)
        if not os.path.exists(p):
            return None
        try:
            with open(p) as fh:
                return json.load(fh)
        except (OSError, json.JSONDecodeError):
            return None

    def put(self, url: str, payload, status: str) -> None:
        if not self.dir:
            return
        try:
            with open(self._path(url), "w") as fh:
                json.dump({"url": url, "status": status, "fetched": _now(),
                           "payload": payload}, fh)
        except OSError:
            pass


# --------------------------------------------------------------------- parsing
@dataclass
class Annotation:
    target: str
    name: str = ""
    pfam: str = ""
    clan: str = ""
    organism: str = ""
    taxon_id: str = ""
    length: str = ""
    status: str = "not_fetched"
    source: str = ""

    def row(self) -> list[str]:
        return [getattr(self, c) for c in ANNOTATION_COLUMNS]


def parse_uniprot(acc: str, payload: dict) -> Annotation:
    """Pull name, Pfam accessions, organism and length from a UniProtKB record.

    Written against the documented UniProtKB JSON shape and deliberately
    forgiving: any field that is missing or differently nested yields an empty
    string rather than an exception, so a shape change degrades the row instead
    of ending the run. Confirm the shape with ``--probe`` before trusting a
    bulk run.
    """
    a = Annotation(target=acc, status="ok", source="uniprotkb")
    if not isinstance(payload, dict):
        return Annotation(target=acc, status="parse_error", source="uniprotkb")

    desc = payload.get("proteinDescription") or {}
    for key in ("recommendedName", "submissionNames"):
        node = desc.get(key)
        node = node[0] if isinstance(node, list) and node else node
        if isinstance(node, dict):
            full = (node.get("fullName") or {}).get("value")
            if full:
                a.name = full
                break
    if not a.name:
        a.name = payload.get("uniProtkbId") or ""

    org = payload.get("organism") or {}
    a.organism = org.get("scientificName") or ""
    tax = org.get("taxonId")
    a.taxon_id = str(tax) if tax is not None else ""

    seq = payload.get("sequence") or {}
    ln = seq.get("length")
    a.length = str(ln) if ln is not None else ""

    pfams = []
    for xref in payload.get("uniProtKBCrossReferences") or []:
        if isinstance(xref, dict) and xref.get("database") == "Pfam":
            pid = xref.get("id")
            if pid:
                pfams.append(pid)
    a.pfam = ";".join(dict.fromkeys(pfams))

    if payload.get("entryType", "").lower().startswith("inactive") or payload.get("inactiveReason"):
        a.status = "withdrawn"
    return a


def parse_uniparc(acc: str, payload: dict) -> Annotation:
    """Recover organism and protein name for a record no longer in UniProtKB."""
    a = Annotation(target=acc, status="withdrawn_recovered", source="uniparc")
    results = payload.get("results") if isinstance(payload, dict) else None
    if not results:
        return Annotation(target=acc, status="not_found", source="uniparc")
    entry = results[0] if isinstance(results, list) else results
    seq = (entry or {}).get("sequence") or {}
    ln = seq.get("length")
    a.length = str(ln) if ln is not None else ""
    for xref in (entry or {}).get("uniParcCrossReferences") or []:
        if not isinstance(xref, dict):
            continue
        if not a.name:
            a.name = xref.get("proteinName") or ""
        org = xref.get("organism") or {}
        if not a.organism:
            a.organism = org.get("scientificName") or ""
            tax = org.get("taxonId")
            a.taxon_id = str(tax) if tax is not None else ""
        if a.name and a.organism:
            break
    return a


def parse_interpro_clan(payload: dict) -> str:
    """Clan (set) accession for a Pfam family, or an empty string."""
    if not isinstance(payload, dict):
        return ""
    meta = (payload.get("metadata") or {})
    sets = meta.get("set_info") or meta.get("sets") or {}
    if isinstance(sets, dict):
        return sets.get("accession") or ""
    if isinstance(sets, list) and sets:
        first = sets[0]
        if isinstance(first, dict):
            return first.get("accession") or ""
        if isinstance(first, str):
            return first
    return ""


# ------------------------------------------------------------------- pipeline
@dataclass
class Manifest:
    created: str = field(default_factory=_now)
    endpoints: list[str] = field(default_factory=list)
    requested: int = 0
    counts: dict = field(default_factory=dict)
    output: str = ""
    output_sha256: str = ""
    note: str = ("Every requested accession has a row. A non-ok status is a recorded "
                 "retrieval outcome, never a biological absence.")


def annotate(accessions: list[str], cache_dir: str | None = None,
             clans: bool = True, pause: float = 0.34,
             fetcher=fetch_json) -> tuple[list[Annotation], Manifest]:
    """Build one annotation row per accession. ``fetcher`` is injectable for tests."""
    if not accessions:
        raise AnnotateError("no accessions given")
    cache = Cache(cache_dir)
    out: list[Annotation] = []
    clan_of: dict[str, str] = {}
    endpoints = {UNIPROT}

    for acc in accessions:
        url = UNIPROT.format(acc=acc)
        hit = cache.get(url)
        if hit is not None:
            payload, status = hit.get("payload"), hit.get("status", "ok")
        else:
            payload, status = fetcher(url)
            cache.put(url, payload, status)

        if status == "ok" and payload:
            a = parse_uniprot(acc, payload)
        else:
            a = Annotation(target=acc, status=status, source="uniprotkb")

        if a.status in ("not_found", "withdrawn"):
            endpoints.add(UNIPARC_BY_ACC)
            u2 = UNIPARC_BY_ACC.format(acc=urllib.parse.quote(acc))
            hit2 = cache.get(u2)
            if hit2 is not None:
                p2, s2 = hit2.get("payload"), hit2.get("status", "ok")
            else:
                p2, s2 = fetcher(u2)
                cache.put(u2, p2, s2)
            if s2 == "ok" and p2:
                recovered = parse_uniparc(acc, p2)
                if recovered.status == "withdrawn_recovered":
                    recovered.pfam = a.pfam
                    a = recovered

        if clans and a.pfam:
            for pf in a.pfam.split(";"):
                if pf and pf not in clan_of:
                    endpoints.add(INTERPRO_PFAM)
                    u3 = INTERPRO_PFAM.format(pfam=pf)
                    hit3 = cache.get(u3)
                    if hit3 is not None:
                        p3, s3 = hit3.get("payload"), hit3.get("status", "ok")
                    else:
                        p3, s3 = fetcher(u3)
                        cache.put(u3, p3, s3)
                    clan_of[pf] = parse_interpro_clan(p3) if s3 == "ok" and p3 else ""
            a.clan = ";".join(dict.fromkeys(
                c for c in (clan_of.get(pf, "") for pf in a.pfam.split(";")) if c))
        out.append(a)

    counts: dict[str, int] = {}
    for a in out:
        counts[a.status] = counts.get(a.status, 0) + 1
    man = Manifest(endpoints=sorted(endpoints), requested=len(accessions), counts=counts)
    return out, man


def write_annotations(rows: list[Annotation], path: str, delimiter: str = "\t") -> str:
    with open(path, "w") as fh:
        fh.write(delimiter.join(ANNOTATION_COLUMNS) + "\n")
        for a in rows:
            fh.write(delimiter.join(str(v).replace(delimiter, " ") for v in a.row()) + "\n")
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(man: Manifest, path: str) -> None:
    with open(path, "w") as fh:
        json.dump(asdict(man), fh, indent=2)
        fh.write("\n")


def read_accessions(path: str) -> list[str]:
    """One accession per line, or the first column of a TSV."""
    out: list[str] = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            tok = line.strip().split("\t")[0].strip()
            if not tok or tok.startswith("#"):
                continue
            if i == 0 and tok.lower() in ("target", "candidate", "candidate_id", "accession"):
                continue
            out.append(tok)
    if not out:
        raise AnnotateError(f"no accessions in {path}")
    return list(dict.fromkeys(out))
