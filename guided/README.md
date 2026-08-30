# FOLDTRACE Guided

A paced, checkpointed wrapper around FOLDTRACE for undergraduate research. One
project studies **one candidate structure** against an experimentally characterized
reference and its prespecified chemistry-defining sites, moving through five stages:

| stage | you do | the tool does |
|---|---|---|
| **observe** | read the reference chemistry and the fold-search evidence | prints a briefing of the sites + search metrics |
| **predict** | commit, in writing, to retained / altered / lost per site | records and validates your hypothesis |
| **compute** | run the mapping | superposes the candidate, reads each site, **locks your prediction** |
| **compare** | see where you were right and where you were surprised | scores prediction vs result; unresolved sites are not counted against you |
| **decide** | write a conclusion and the next experiment | records it |

Two rules make it a teaching tool, not just a script:

- **Checkpoints** — a stage will not run until the one before it is done.
- **Prediction lock** — once you `compute`, your prediction is frozen. You cannot
  rewrite the hypothesis after seeing the answer, which is what makes `compare` honest.

The science is unchanged: Guided calls the released `foldtrace.mapping.map_candidate`
with the same TM gate (0.5) and CA-offset threshold (4.0 A), and the same
retained / altered / lost / unresolved definitions. FOLDTRACE tracks *prespecified*
positions from a known reference; it does not discover catalytic residues de novo.

## Command line

A project is a single JSON journal, so students can stop and resume across sessions.

```bash
foldtrace guided init    --project my.json --name my-first-tir \
    --reference examples/tir/reference/SARM1_TIR_6O0R_A.pdb \
    --sites     examples/tir/catalytic_sites.tsv \
    --candidate examples/tir/candidates/A0A933QTC5.pdb
foldtrace guided observe --project my.json --notes "Glu642 is the catalytic key"
foldtrace guided predict --project my.json \
    --set Tyr568=retained --set Trp638=retained --set Glu642=lost \
    --rationale "Strong fold match but the catalytic Glu looks gone"
foldtrace guided compute --project my.json    # locks the prediction
foldtrace guided compare --project my.json
foldtrace guided decide  --project my.json --conclusion "Non-functional; Glu642 lost" \
    --next "Screen A0A953THP3 (predicted retained)"
foldtrace guided report  --project my.json    # Markdown write-up
```

`guided/examples/run_guided_demo.sh` runs exactly this against the bundled TIR example.

## Python

```python
from foldtrace.guided import GuidedProject
p = GuidedProject("my-first-tir", ref_pdb, sites_tsv, candidate_pdb)
p.observe()
p.predict({"Tyr568": "retained", "Trp638": "retained", "Glu642": "lost"})
result = p.compute()          # returns the released CandidateResult; locks the prediction
scorecard = p.compare()       # {'n_correct': ..., 'surprises': [...], ...}
p.decide("Non-functional; Glu642 lost")
print(p.report())
p.save("my.json")             # resume later with GuidedProject.load("my.json")
```
