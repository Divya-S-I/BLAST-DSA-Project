# BLAST Lab: Interactive Sequence Alignment Explorer

A four-page academic Streamlit site:

1. **How BLAST Works**: interactive exact DNA k-mers, query lookup, seed locations, ungapped X-drop extension and ranking.
2. **Redundant Seed Extensions**: preserves the original `redundancy_engine.py` algorithm and original experiment controls, seed tables, interval charts, k sweep and CSV downloads. Style only is adjusted.
3. **BLAST vs Smith-Waterman**: analyse the **192 supplied protein comparison outputs** by family and pair; optional live *short DNA* example comparing heuristic ungapped extension with exact Smith-Waterman.

4. **Adaptive k-mer Selection**: probability DP, sensitivity target, chosen k, random-hit estimates and reproducible simulated trials. Includes the original NCBI BLASTN scripts as supplementary material.

## Run on Mac

```bash
cd /path/to/BLAST_Explorer
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 -m streamlit run How_BLAST_Works.py
```

Open `http://localhost:8501`. Stop with Control+C. Later, repeat only `cd`, `source .venv/bin/activate`, `python3 -m streamlit run How_BLAST_Works.py`.

## Deploy

Create a **new** GitHub repository and upload this folder's contents (excluding `.venv`). Deploy with Streamlit Community Cloud, entrypoint `How_BLAST_Works.py`. This is independent of the earlier BLAST site.

## Notes on validity

- The redundancy implementation is an **educational ungapped exact-word** heuristic, not NCBI BLAST. Coverage-based skipping is **not proven lossless**; score disagreements can occur.
- Very short k on large genome can generate huge hit lists. The experiment caps collected hits; truncated runs are not full coverage.
- The E. coli dataset contains exact 500 bp fragments from its own reference and is a **positive control**, not a realistic sensitivity study.
- The protein results were **precomputed** by a separately supplied notebook; the chart page reads those CSVs rather than rerunning BLASTP. Check scoring consistency, because some BLAST-derived raw scores exceed Smith-Waterman raw scores.
- Modules 1 and 2 use DNA A/C/G/T; Module 3 saved results concern **protein** sequences. Do not equate their scoring systems.
- Reference locations are 0-based/end-exclusive in the educational DNA code; some protein result source columns are 1-based.

## Notes on adaptive-k source scripts
The two bundled scripts require NCBI BLAST+ executables, an NCBI email and access to NCBI data. The Streamlit adaptive page itself runs locally using the mathematical model, **not** those expensive BLASTN experiments; no experimental BLASTN conclusions are pre-populated. The original scripts are unmodified and are under `research_scripts/`.

## Real five-step adaptive BLASTN experiment (Page 4)

The adaptive-k page now has **two distinct modes**:

1. The existing instant mathematical DP model and Bernoulli simulation (no BLAST+ needed).
2. **Actual five-step BLASTN experiment**: generate planted queries from a bundled, pasted, or uploaded single-record DNA reference; run a k=7 BLASTN pilot on five queries; estimate identity from best hit per pilot query; calculate exact Bernoulli k-run probabilities; select the largest k in 7–15 predicted to meet 90% sensitivity; search all queries with that k; compare with fixed k=7,11,15 and a seeded random k, showing measured planted-locus recovery, runtime and CSV export.

This live experiment uses the NCBI **BLAST+ executables**, which are separate from the Python package `biopython`. On macOS with Homebrew installed, run:

```bash
brew install blast
blastn -version
makeblastdb -version
```

Then activate the site's Python virtual environment and restart Streamlit as above. Without BLAST+, Page 4's mathematical model remains available and the live-run button is disabled. The bundled E. coli genome is available offline; the live experiment makes no NCBI network requests.

**Interpretation:** This is a controlled *substitution-only* planted-locus experiment, not general biological sensitivity. A pilot detects only some matches and may bias estimated identity. If no pilot hits are found, the script labels its p=0.70 fallback. The random baseline is a single seeded choice. Timings are one-run estimates; adaptive k includes the pilot time but all methods exclude the shared database construction time. Repeat for robust performance claims.
