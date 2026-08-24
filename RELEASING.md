# Releasing foldtrace and minting a Zenodo DOI

The DOI is created by Zenodo from a GitHub release; it cannot be minted from the repo alone.

1. Sign in at https://zenodo.org with GitHub, open **Settings -> GitHub**, and flip the
   `headsparkle/foldtrace` switch to **On**. (`.zenodo.json` in this repo supplies the metadata.)
2. On GitHub, **Releases -> Draft a new release**, tag `v0.1.0`, title `foldtrace 0.1.0`, publish.
3. Zenodo archives the release automatically and issues a version DOI plus a concept DOI
   (the concept DOI always resolves to the latest version). Copy the concept DOI.
4. Add the DOI badge to the top of `README.md` and the DOI to `CITATION.cff`
   (`identifiers:` / `doi:`), then reference it in the manuscript's Data & Code Availability.

Cite the concept DOI in the paper so the link stays valid across future versions.
