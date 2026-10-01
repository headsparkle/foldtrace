"""
Multi-query Foldseek search for AChE-fold expansion.
Submits 4 bacterial starting structures to Foldseek web API concurrently,
pools results, deduplicates by UniProt accession, writes merged raw TSV.

Queries:
  1QE3   Bacillus subtilis PNBE         Firmicutes / Bacillales
  7X8L   Altericroceibacterium indicum  Alphaproteobacteria
  8S9J   Staphylococcus aureus FphA     Firmicutes / Staphylococcales
  Q73QC2 Treponema denticola (AF model) Spirochaetes — Cys-nucleophile query

Usage:
  python multi_query_foldseek.py
"""
import json, time, sys, gzip, threading
import urllib.request

API       = "https://search.foldseek.com/api"
MAX_HITS  = 3000
OUT_DIR   = "/Users/lukebegg/PanD/ache_project/results"
POOL_TSV  = f"{OUT_DIR}/pooled_raw.tsv"

QUERIES = [
    ("1QE3",   "/Users/lukebegg/PanD/ache_project/reference/1QE3.pdb"),
    ("7X8L",   "/Users/lukebegg/PanD/ache_project/reference/7X8L.pdb"),
    ("8S9J",   "/Users/lukebegg/PanD/ache_project/reference/8S9J.pdb"),
    ("Q73QC2", "/Users/lukebegg/PanD/ache_project/models/Q73QC2.pdb"),
]

def submit_job(name, pdb_path):
    with open(pdb_path, "rb") as f:
        pdb_bytes = f.read()
    fname = pdb_path.split("/")[-1]
    boundary = "----FormBoundary7MA4YWxkTrZu0gW"
    body = (
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="q"; filename="{fname}"\r\n'
        f"Content-Type: application/octet-stream\r\n\r\n"
    ).encode() + pdb_bytes + (
        f"\r\n--{boundary}\r\n"
        f'Content-Disposition: form-data; name="mode"\r\n\r\n'
        f"3diaa\r\n"
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="database[]"\r\n\r\n'
        f"afdb50\r\n"
        f"--{boundary}--\r\n"
    ).encode()
    req = urllib.request.Request(
        f"{API}/ticket",
        data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        ticket = json.loads(r.read())
    print(f"  [{name}] submitted → ticket {ticket['id']}", flush=True)
    return ticket["id"]

def poll_and_download(name, ticket_id, results):
    out_path = f"{OUT_DIR}/raw_{name}.tsv"
    for attempt in range(180):  # max 15 min
        time.sleep(5)
        with urllib.request.urlopen(f"{API}/ticket/{ticket_id}", timeout=30) as r:
            status = json.loads(r.read())
        state = status.get("status", "")
        if attempt % 6 == 0:
            print(f"  [{name}] {attempt*5}s — {state}", flush=True)
        if state == "COMPLETE":
            break
        if state == "ERROR":
            print(f"  [{name}] ERROR: {status}")
            results[name] = None
            return
    else:
        print(f"  [{name}] TIMED OUT")
        results[name] = None
        return

    with urllib.request.urlopen(f"{API}/result/download/{ticket_id}", timeout=120) as r:
        data = r.read()
    try:
        raw = gzip.decompress(data).decode("utf-8")
    except (gzip.BadGzipFile, OSError):
        raw = data.decode("utf-8")

    lines = [l for l in raw.splitlines() if l.strip()]
    with open(out_path, "w") as f:
        f.write(raw)
    print(f"  [{name}] DONE — {len(lines)} hits → {out_path}", flush=True)
    results[name] = (out_path, lines)


def extract_accession(target_field):
    """Extract UniProt accession from 'AF-A0A123-F1-model_v6 ...' format."""
    tok = target_field.split()[0]          # AF-A0A123-F1-model_v6
    parts = tok.split("-")
    if len(parts) >= 3 and parts[0] == "AF":
        return parts[1]
    return tok                             # fallback: return as-is

# ── Step 1: submit all 4 jobs ─────────────────────────────────────────
print("=== Submitting 4 Foldseek jobs ===", flush=True)
tickets = {}
for name, path in QUERIES:
    try:
        tid = submit_job(name, path)
        tickets[name] = tid
    except Exception as e:
        print(f"  [{name}] submit failed: {e}")

# ── Step 2: poll and download in parallel threads ─────────────────────
print("\n=== Polling for results (up to 15 min each) ===", flush=True)
results = {}
threads = []
for name, tid in tickets.items():
    t = threading.Thread(target=poll_and_download, args=(name, tid, results))
    t.start()
    threads.append(t)
for t in threads:
    t.join()

# ── Step 3: pool and deduplicate ──────────────────────────────────────
print("\n=== Pooling and deduplicating ===", flush=True)
seen_accessions = {}          # accession → (name, line, score)
all_lines_by_query = {}

for name, val in results.items():
    if val is None:
        print(f"  Skipping {name} (failed)")
        continue
    out_path, lines = val
    all_lines_by_query[name] = lines
    for line in lines:
        cols = line.split("\t")
        if len(cols) < 2:
            continue
        acc = extract_accession(cols[1])
        # Keep the hit from any query; prefer higher bit-score (col 11 in m8)
        try:
            bits = float(cols[11]) if len(cols) > 11 else 0.0
        except ValueError:
            bits = 0.0
        if acc not in seen_accessions or bits > seen_accessions[acc][2]:
            seen_accessions[acc] = (name, line, bits)

print(f"  Total unique accessions: {len(seen_accessions)}", flush=True)

# ── Step 4: write pooled TSV ──────────────────────────────────────────
pooled_lines = [row[1] for row in seen_accessions.values()]
with open(POOL_TSV, "w") as f:
    f.write("\n".join(pooled_lines) + "\n")
print(f"  Pooled TSV → {POOL_TSV}  ({len(pooled_lines)} lines)", flush=True)

# ── Step 5: summary ───────────────────────────────────────────────────
print("\n=== Per-query hit counts ===")
for name, val in results.items():
    if val:
        print(f"  {name}: {len(val[1])} raw hits")
    else:
        print(f"  {name}: FAILED")
print(f"\n  After dedup: {len(seen_accessions)} unique accessions")
print(f"\nNext step:")
print(f"  source /Users/lukebegg/PanD/foldtrace/.venv/bin/activate")
print(f"  foldtrace search --from-foldseek {POOL_TSV} --max-hits {MAX_HITS} --out {OUT_DIR}/pooled_hits.tsv")
