"""The GuidedProject research notebook that paces a FOLDTRACE investigation.

FOLDTRACE Guided turns a single structure-first active-site investigation into an
ordered, checkpointed lab notebook aimed at undergraduate researchers. One project
studies one candidate structure against one experimentally characterized reference
and its prespecified chemistry-defining sites, moving through five stages:

    observe  -> look at the reference chemistry and the fold-search evidence
    predict  -> commit, in writing, to what each site will be (retained/altered/lost)
    compute  -> run the released FOLDTRACE mapping (this LOCKS the prediction)
    compare  -> score the prediction against what FOLDTRACE actually read
    decide   -> record a conclusion and the next experiment

The pedagogy lives in two rules the class enforces:

* **Checkpoints.** A stage cannot run before the stage before it is done
  (:class:`CheckpointError`). The workflow is the point, not just the answer.
* **Prediction lock.** Once you compute, your prediction is frozen
  (:class:`PredictionLockedError`); you cannot quietly rewrite the hypothesis after
  seeing the result. This is what makes ``compare`` an honest self-assessment.

The scientific behaviour, thresholds, and state definitions are exactly those of the
released tool: this module calls :func:`foldtrace.mapping.map_candidate` and changes
nothing about how a site is called retained / altered / lost / unresolved. FOLDTRACE
tracks *prespecified* positions from a known reference; it does not discover catalytic
residues de novo, and Guided does not either.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone

from ..io import load_structure, parse_sites
from ..mapping import (
    map_candidate, CandidateResult, SiteCall,
    RETAINED, ALTERED, LOST, UNRESOLVED,
)

# Thresholds re-exported so a guided project shares the released defaults verbatim.
TM_GATE = 0.5      # foldtrace.mapping.map_candidate default tm_gate
OFFSET_A = 4.0     # foldtrace.mapping.map_candidate default offset_threshold

# The states a student may predict (unresolved is an outcome of the method, not a
# hypothesis a student commits to, so it is not offered as a prediction choice).
PREDICTABLE_STATES = (RETAINED, ALTERED, LOST)


@dataclass(frozen=True)
class Stage:
    id: str
    title: str
    prompt: str      # the guiding question shown to the student at this stage


STAGES: tuple[Stage, ...] = (
    Stage("observe", "Observe",
          "What is the reference's chemistry, and how good is this candidate's fold match? "
          "Read the prespecified sites and the fold-search evidence before forming any opinion."),
    Stage("predict", "Predict",
          "For each site, commit in writing to retained, altered, or lost, with your reasoning. "
          "You are locking in a hypothesis before you look at the structural answer."),
    Stage("compute", "Compute",
          "Run FOLDTRACE mapping. This superposes the candidate on the reference and reads each "
          "site's state. Running it freezes your prediction."),
    Stage("compare", "Compare",
          "Where did FOLDTRACE agree with you, and where were you surprised? Every surprise is a "
          "place to look harder at the structure."),
    Stage("decide", "Decide",
          "Given the evidence, what do you conclude about this candidate, and what is the next "
          "experiment (a wet-lab assay, another candidate, a tighter threshold)?"),
)
STAGE_IDS: tuple[str, ...] = tuple(s.id for s in STAGES)
_STAGE_INDEX = {s.id: i for i, s in enumerate(STAGES)}


class GuidedError(Exception):
    """Base class for guided-workflow errors."""


class CheckpointError(GuidedError):
    """Raised when a stage is attempted before the stage before it is complete."""


class PredictionLockedError(GuidedError):
    """Raised when a prediction is changed after ``compute`` has locked it."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _sitecall_to_dict(c: SiteCall) -> dict:
    return {"label": c.label, "ref_resnum": c.ref_resnum, "obs_resnum": c.obs_resnum,
            "obs_res": c.obs_res, "ca_offset": c.ca_offset, "state": c.state, "reason": c.reason}


def _result_to_dict(r: CandidateResult) -> dict:
    return {"name": r.name, "tm_norm_ref": r.tm_norm_ref, "rmsd": r.rmsd, "fold_ok": r.fold_ok,
            "verdict": r.verdict, "state_reason": r.state_reason,
            "sites": [_sitecall_to_dict(c) for c in r.sites]}


class GuidedProject:
    """One paced FOLDTRACE investigation of a single candidate, persisted as a JSON journal.

    Parameters
    ----------
    name : a short label for the investigation.
    reference_path, sites_path, candidate_path : the released-tool inputs (a reference
        structure, its prespecified-sites TSV, and the one candidate structure to study).
    hit : optional dict of fold-search metrics for the candidate (e.g. ``prob``,
        ``fident``, ``alntmscore``) surfaced during ``observe`` so predictions rest on
        the search evidence.
    tm_gate, offset_threshold : passed straight through to
        :func:`foldtrace.mapping.map_candidate`; default to the released values.
    """

    def __init__(self, name: str, reference_path: str, sites_path: str, candidate_path: str,
                 *, hit: dict | None = None, tm_gate: float = TM_GATE,
                 offset_threshold: float = OFFSET_A):
        self.name = name
        self.reference_path = reference_path
        self.sites_path = sites_path
        self.candidate_path = candidate_path
        self.hit = dict(hit) if hit else None
        self.tm_gate = tm_gate
        self.offset_threshold = offset_threshold
        self.created = _now()
        self.sites = parse_sites(sites_path)   # validates inputs early; exposes site labels
        self.stages: dict[str, dict] = {sid: {"done": False} for sid in STAGE_IDS}

    # -- helpers ----------------------------------------------------------------
    @property
    def site_labels(self) -> list[str]:
        return [s.label for s in self.sites]

    def _require(self, stage_id: str) -> None:
        """Raise CheckpointError unless every stage before ``stage_id`` is done."""
        idx = _STAGE_INDEX[stage_id]
        for prior in STAGE_IDS[:idx]:
            if not self.stages[prior]["done"]:
                nice = STAGES[_STAGE_INDEX[prior]].title
                raise CheckpointError(
                    f"cannot start '{stage_id}': the '{prior}' ({nice}) stage is not done yet")

    def _mark_done(self, stage_id: str, **payload) -> None:
        self.stages[stage_id].update(done=True, timestamp=_now(), **payload)

    # -- stage 1: observe -------------------------------------------------------
    def briefing(self) -> dict:
        """The reference chemistry + fold-search evidence a student weighs before predicting."""
        return {
            "reference": os.path.basename(self.reference_path),
            "candidate": os.path.splitext(os.path.basename(self.candidate_path))[0],
            "tm_gate": self.tm_gate,
            "offset_threshold": self.offset_threshold,
            "fold_search": self.hit,
            "sites": [{"label": s.label, "ref_resnum": s.resnum, "role": s.role,
                       "expected": sorted(s.expected), "altered": sorted(s.altered),
                       "key": s.key} for s in self.sites],
        }

    def observe(self, notes: str = "") -> dict:
        """Record what the student sees. Returns (and stores) the briefing."""
        brief = self.briefing()
        self._mark_done("observe", notes=notes, briefing=brief)
        return brief

    # -- stage 2: predict -------------------------------------------------------
    def predict(self, site_predictions: dict[str, str], rationale: str = "",
                verdict: str | None = None) -> dict:
        """Commit a predicted state (retained/altered/lost) for each site.

        May be revised freely until :meth:`compute` locks it, after which it raises
        :class:`PredictionLockedError`.
        """
        self._require("predict")
        if self.stages["predict"].get("locked"):
            raise PredictionLockedError(
                "prediction is locked (compute has run); a hypothesis cannot change after the result")

        labels = set(self.site_labels)
        given = set(site_predictions)
        if given - labels:
            raise ValueError(f"unknown site labels in prediction: {sorted(given - labels)}")
        if labels - given:
            raise ValueError(f"missing predictions for sites: {sorted(labels - given)}")
        bad = {lab: st for lab, st in site_predictions.items() if st not in PREDICTABLE_STATES}
        if bad:
            raise ValueError(
                f"predictions must be one of {PREDICTABLE_STATES}; got {bad}")
        if verdict is not None and verdict not in PREDICTABLE_STATES:
            raise ValueError(f"verdict must be one of {PREDICTABLE_STATES} or None; got {verdict!r}")

        record = {"site_predictions": dict(site_predictions),
                  "verdict": verdict, "rationale": rationale}
        self._mark_done("predict", locked=False, **record)
        return record

    # -- stage 3: compute -------------------------------------------------------
    def compute(self) -> CandidateResult:
        """Run the released FOLDTRACE mapping and lock the prediction."""
        self._require("compute")
        ref = load_structure(self.reference_path)
        cand = load_structure(self.candidate_path)
        candidate_name = os.path.splitext(os.path.basename(self.candidate_path))[0]
        result = map_candidate(ref, cand, self.sites, candidate_name=candidate_name,
                               offset_threshold=self.offset_threshold, tm_gate=self.tm_gate)
        self.stages["predict"]["locked"] = True
        self._mark_done("compute", result=_result_to_dict(result))
        return result

    # -- stage 4: compare -------------------------------------------------------
    def compare(self) -> dict:
        """Score the locked prediction against the computed states, site by site."""
        self._require("compare")
        pred = self.stages["predict"]["site_predictions"]
        computed = {c["label"]: c["state"] for c in self.stages["compute"]["result"]["sites"]}
        rows = []
        n_correct = n_scored = 0
        for label in self.site_labels:
            predicted = pred[label]
            observed = computed.get(label, UNRESOLVED)
            # unresolved is not scored against the student: the method could not read it
            scorable = observed != UNRESOLVED
            hit = scorable and predicted == observed
            if scorable:
                n_scored += 1
                n_correct += int(hit)
            rows.append({"label": label, "predicted": predicted, "observed": observed,
                         "match": hit, "scored": scorable,
                         "surprise": scorable and not hit})
        result = self.stages["compute"]["result"]
        scorecard = {
            "computed_verdict": result["verdict"],
            "predicted_verdict": self.stages["predict"].get("verdict"),
            "n_scored": n_scored, "n_correct": n_correct,
            "n_unresolved": sum(1 for r in rows if not r["scored"]),
            "surprises": [r["label"] for r in rows if r["surprise"]],
            "sites": rows,
        }
        self._mark_done("compare", scorecard=scorecard)
        return scorecard

    # -- stage 5: decide --------------------------------------------------------
    def decide(self, conclusion: str, next_action: str = "") -> dict:
        """Record the student's conclusion and the next experiment."""
        self._require("decide")
        record = {"conclusion": conclusion, "next_action": next_action}
        self._mark_done("decide", **record)
        return record

    # -- status / reporting -----------------------------------------------------
    def status(self) -> dict:
        done = [sid for sid in STAGE_IDS if self.stages[sid]["done"]]
        nxt = next((sid for sid in STAGE_IDS if not self.stages[sid]["done"]), None)
        return {"name": self.name, "done": done, "next": nxt,
                "complete": nxt is None,
                "prediction_locked": bool(self.stages["predict"].get("locked"))}

    def report(self) -> str:
        """A short Markdown write-up of the investigation so far."""
        st = self.status()
        lines = [f"# FOLDTRACE Guided: {self.name}", ""]
        lines.append(f"Reference `{os.path.basename(self.reference_path)}` vs candidate "
                     f"`{os.path.splitext(os.path.basename(self.candidate_path))[0]}`  ")
        lines.append(f"Stages done: {', '.join(st['done']) or 'none'}"
                     + (f" | next: {st['next']}" if st["next"] else " | complete"))
        if self.stages["predict"]["done"]:
            lines += ["", "## Prediction", ""]
            for lab, state in self.stages["predict"]["site_predictions"].items():
                lines.append(f"- {lab}: **{state}**")
            if self.stages["predict"].get("rationale"):
                lines += ["", f"> {self.stages['predict']['rationale']}"]
        if self.stages["compute"]["done"]:
            r = self.stages["compute"]["result"]
            lines += ["", "## FOLDTRACE result", "",
                      f"Verdict: **{r['verdict']}** (TM {r['tm_norm_ref']:.2f}, RMSD {r['rmsd']:.2f} A)"]
            for c in r["sites"]:
                obs = f"{c['obs_res']}{c['obs_resnum']}" if c["obs_resnum"] is not None else "-"
                lines.append(f"- {c['label']}: {c['state']} ({obs})")
        if self.stages["compare"]["done"]:
            sc = self.stages["compare"]["scorecard"]
            lines += ["", "## Compare", "",
                      f"Correct on {sc['n_correct']}/{sc['n_scored']} scored sites"
                      + (f"; {sc['n_unresolved']} unresolved (not scored)" if sc["n_unresolved"] else "")]
            if sc["surprises"]:
                lines.append(f"Surprises to investigate: {', '.join(sc['surprises'])}")
        if self.stages["decide"]["done"]:
            d = self.stages["decide"]
            lines += ["", "## Decision", "", d["conclusion"]]
            if d.get("next_action"):
                lines += ["", f"**Next:** {d['next_action']}"]
        return "\n".join(lines) + "\n"

    # -- persistence ------------------------------------------------------------
    def to_dict(self) -> dict:
        return {"schema": "foldtrace-guided/1", "name": self.name, "created": self.created,
                "reference_path": self.reference_path, "sites_path": self.sites_path,
                "candidate_path": self.candidate_path, "hit": self.hit,
                "tm_gate": self.tm_gate, "offset_threshold": self.offset_threshold,
                "stages": self.stages}

    def save(self, path: str) -> str:
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=2)
        return path

    @classmethod
    def load(cls, path: str) -> "GuidedProject":
        with open(path) as fh:
            d = json.load(fh)
        proj = cls(d["name"], d["reference_path"], d["sites_path"], d["candidate_path"],
                   hit=d.get("hit"), tm_gate=d.get("tm_gate", TM_GATE),
                   offset_threshold=d.get("offset_threshold", OFFSET_A))
        proj.created = d.get("created", proj.created)
        proj.stages = d["stages"]
        for sid in STAGE_IDS:                 # tolerate journals written by older code
            proj.stages.setdefault(sid, {"done": False})
        return proj
