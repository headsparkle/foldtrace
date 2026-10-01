"""Submit 1EVE (AChE) to Foldseek web API, wait for results, save TSV."""
import json, time, sys, gzip
import urllib.request, urllib.parse

API = "https://search.foldseek.com/api"
PDB_PATH = "/Users/lukebegg/PanD/ache_project/reference/1EVE.pdb"
OUT_PATH  = "/Users/lukebegg/PanD/ache_project/results/foldseek_raw.tsv"

# ── Submit job ────────────────────────────────────────────────────────
print("Reading 1EVE.pdb ...", flush=True)
with open(PDB_PATH, "rb") as f:
    pdb_bytes = f.read()

# Multipart form-data submission
boundary = "----FormBoundary7MA4YWxkTrZu0gW"
body = (
    f"--{boundary}\r\n"
    f'Content-Disposition: form-data; name="q"; filename="1EVE.pdb"\r\n'
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
print("Submitting to Foldseek web API (afdb50, 3diaa+aa) ...", flush=True)
with urllib.request.urlopen(req, timeout=60) as r:
    ticket = json.loads(r.read())
ticket_id = ticket["id"]
print(f"Ticket: {ticket_id}", flush=True)

# ── Poll until done ───────────────────────────────────────────────────
for attempt in range(120):  # max 10 min
    time.sleep(5)
    with urllib.request.urlopen(f"{API}/ticket/{ticket_id}", timeout=30) as r:
        status = json.loads(r.read())
    state = status.get("status", "")
    print(f"  [{attempt*5}s] status: {state}", flush=True)
    if state == "COMPLETE":
        break
    if state == "ERROR":
        print("ERROR:", status)
        sys.exit(1)
else:
    print("Timed out waiting for Foldseek job.")
    sys.exit(1)

# ── Download result ───────────────────────────────────────────────────
print("Downloading results ...", flush=True)
with urllib.request.urlopen(f"{API}/result/download/{ticket_id}", timeout=60) as r:
    data = r.read()
    try:
        raw = gzip.decompress(data).decode("utf-8")
    except (gzip.BadGzipFile, OSError):
        raw = data.decode("utf-8")

with open(OUT_PATH, "w") as f:
    f.write(raw)

lines = [l for l in raw.splitlines() if l.strip()]
print(f"Done — {len(lines)} hits written to {OUT_PATH}")
