# Releasing foldtrace

The Zenodo **concept DOI is already minted: `10.5281/zenodo.22081699`** (it always resolves to
the latest version). It is recorded in `README.md`, `CITATION.cff`, and `.zenodo.json`, and is the
DOI cited in the manuscript. Each GitHub release adds a new *version* DOI under that concept DOI.

## v0.2.0 — do NOT release yet

The v0.2.0 changes are staged and tested but must not be released until the author approves.

When approved:

1. Confirm the Zenodo <-> GitHub integration is **On** for `headsparkle/foldtrace`
   (https://zenodo.org, Settings -> GitHub).
2. On GitHub, **Releases -> Draft a new release**, tag `v0.2.0`, title `foldtrace 0.2.0`, publish.
3. Zenodo archives the release automatically and issues the v0.2.0 version DOI under the
   existing concept DOI `10.5281/zenodo.22081699`.

No manuscript edit is needed on release: the paper cites the concept DOI, which is stable across
versions.
