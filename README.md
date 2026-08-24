# foldsite

**Structure-first active-site mapping: read the catalytic state of a fold from predicted structures.**

Sequence homology assigns a protein the function of its closest annotated relative. That
fails when a familiar fold is retained after the chemistry inside it has been altered,
repurposed, or lost. `foldsite` reads the active site directly from structure instead: it
superposes a candidate on a reference with TM-align and reports, for each chemistry-defining
position, whether the catalytic residue is **retained**, **lost**, or **unresolved**.

The reading is *order-aware*. Rather than taking the residue whose Cα happens to sit nearest a
reference catalytic atom (which can grab a spatially-close but sequence-unrelated residue),
`foldsite` uses the TM-align residue correspondence, which preserves sequence order, and only
then checks that the corresponding residue lies within a distance cutoff of the reference
position. A call is emitted only when fold correspondence and spatial placement agree.

This is the active-site-mapping stage of the workflow in Zimmer, *Structure-first search and
active-site mapping across protein folds* (see `CITATION.cff`). v0.1 implements that stage as a
standalone, installable tool with a worked TIR demo; the upstream fold **search** stage is a
thin wrapper over [Foldseek](https://github.com/steineggerlab/foldseek) (see below).

## Install

```bash
git clone https://github.com/headsparkle/foldsite
cd foldsite
pip install -e .
```

Dependencies are pinned in `requirements.txt` (NumPy, Biopython 1.82, tmtools 0.3.0). TM-align
runs in-process via `tmtools`; no external binary is required. Python >= 3.9.

## Worked demo (TIR NADase, ~1 s)

```bash
foldsite map \
  --reference examples/tir/reference/SARM1_TIR_6O0R_A.pdb \
  --sites     examples/tir/catalytic_sites.tsv \
  --candidates examples/tir/candidates/*.pdb \
  --out       out.tsv
```

or just `bash examples/tir/run_demo.sh`. Reference: the SARM1 TIR domain (PDB 6O0R chain A,
residues 561-700) with its Tyr568 / Trp638 / Glu642 NADase site. Candidates: four AlphaFold
models bundled in `examples/tir/candidates/`. Runtime is about **1 second for the four
candidates on a single CPU core**. Expected output (also checked in as
`examples/tir/expected_output.tsv`):

| candidate | tm_norm_ref | Glu642_obs | Glu642_state | verdict |
|---|---|---|---|---|
| A0A953THP3 | 0.886 | E147 | retained | **retained** |
| A0A7Y7TLD0 | 0.898 | E83  | retained | **retained** |
| Q9FHM1     | 0.801 | E84  | retained | **retained** (a TIR NADase site currently annotated as a CD38-like enzyme) |
| A0A933QTC5 | 0.752 | T80  | lost     | **lost** (catalytic Glu → Thr; pocket Trp → Leu76) |

A0A933QTC5 is the naturally catalytic-Glu-lost member; the others retain the full site. These
reproduce the calls in the manuscript.

## Input

- **`--reference`**: a PDB/mmCIF structure whose catalytic residues you know.
- **`--sites`**: a TSV of the chemistry-defining positions:

  ```
  label    ref_resnum    expected    role               key
  Tyr568   568           Y,F         pocket_aromatic    0
  Trp638   638           W           nicotinamide_stack 0
  Glu642   642           E           catalytic          1
  ```

  `expected` lists the one-letter residues that count as *retained* (several codes express a
  conservative/altered-but-functional set, e.g. `Y,F`). `key=1` marks the residue that drives
  the overall verdict.
- **`--candidates`**: one or more structures (files or globs).

Tuning: `--offset-threshold` (default 4.0 Å) is the maximum Cα-Cα offset for a trusted call;
`--tm-gate` (default 0.5) is the minimum reference-normalised TM-score for the fold to be
treated as the same.

## The search stage

The discovery step that produces candidates is a fold search with Foldseek against the
AlphaFold database, followed by removal of already-annotated family members. v0.1 documents
this rather than re-implementing it; the canonical call is:

```bash
foldseek easy-search reference.pdb afdb50 hits.m8 tmp --format-output "query,target,evalue,alntmscore"
```

Feed the resulting hit structures to `foldsite map`. A native `foldsite search` wrapper is
planned for v0.2.

## Tests

```bash
pip install -e ".[test]"
pytest
```

## Citing

See `CITATION.cff`. TM-align: Zhang & Skolnick, *Nucleic Acids Res.* 33, 2302 (2005);
`tmtools` Python bindings: https://github.com/jvkersch/tmtools; Foldseek: van Kempen et al.,
*Nat. Biotechnol.* 42, 243 (2024).

## License

MIT (see `LICENSE`).
