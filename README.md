

# Crohn's Disease Meta-Analyses — Temporal & Causal Microbiome Pipeline

Overview
--------
This repository contains a reproducible analysis pipeline for temporal and causal meta-analyses of gut microbiome data from Crohn's disease (IBD) cohorts. The project aims to:
- Preprocess and clean microbiome abundance tables from longitudinal studies.
- Perform differential abundance and temporal correlation analyses.
- Apply causal discovery methods (DBN, Granger causality, Tigramite) to infer directed relationships and temporal causality between microbial taxa and clinical states (e.g., flare vs remission).
- Produce consensus causal edges and validation-ready outputs for follow-up analysis and manuscript figures.

Repository structure
--------------------
- data — Raw and processed input CSVs:
  - `hmp2_ibd_metagenomics_atlas_20260219_121629.csv` (original atlas)
  - `hmp2_ibd_metagenomics_cleaned.csv` (cleaned table)
  - `hmp2_ibd_metagenomics_preprocessed.csv` (normalized / preprocessed for analysis)
  - `hmp2_ibd_metagenomics_temporal.csv` (long-format temporal table for time-series analyses)
- src — Analysis scripts and pipeline modules:
  - `clean_dataset.py` — cleaning and QC utilities
  - `preprocess_pipeline.py` — normalization, filtering, and formatting steps
  - `statistical_analysis.py` — differential abundance and summary statistics
  - `temporal_analysis.py` — temporal correlations and lag analysis
  - `causal_discovery_dbn.py` — Dynamic Bayesian Network causal discovery
  - `causal_discovery_granger.py` — Granger causality tests across lags
  - `causal_discovery_tigramite.py` — Tigramite-based PC/PCMCI temporal causal discovery
  - `causal_consensus.py` — consensus aggregation across causal methods
- results — Generated outputs and notes:
  - `differential_abundance_results.csv` — differential taxa tests
  - `temporal_lagged_correlations.csv` — Spearman / lagged correlations
  - `causal_dbn_results.csv`, `causal_granger_results.csv`, `causal_tigramite_results.csv` — method-specific causal outputs
  - `causal_consensus_edges.csv` — aggregated consensus edges across methods
  - `flare_vs_remission_results.csv` — comparisons of microbial signatures across disease states
  - `manuscript_draft.md` and `phase_*` documents — analysis notes and write-ups

Requirements
------------
- Python 3.9+ recommended
- Common Python libraries used (install into a virtualenv):
  - pandas, numpy, scipy, scikit-learn, statsmodels
  - networkx, matplotlib, seaborn
  - tigramite (if using Tigramite analyses)
  - pomegranate or bnlearn (if Dynamic Bayesian Network implementation depends on one)
- Create a virtual environment and install dependencies:
  - python -m venv .venv
  - Windows:
    - .venv\Scripts\Activate.ps1
  - Unix:
    - source .venv/bin/activate
  - pip install -r requirements.txt
Note: If a `requirements.txt` does not exist, generate one with the packages above appropriate to your environment.

How the pipeline is organized (high level)
-----------------------------------------
1. Raw input in data is loaded and QC'd using `clean_dataset.py`.
2. `preprocess_pipeline.py` performs filtering (prevalence/abundance thresholds), normalization (CLR, relative abundance, or log transformation), and converts the dataset to the long temporal format saved as `hmp2_ibd_metagenomics_temporal.csv`.
3. `statistical_analysis.py` runs cross-sectional differential abundance tests (e.g., Wilcoxon/Mann–Whitney, linear models) and writes `differential_abundance_results.csv`.
4. `temporal_analysis.py` computes pairwise lagged correlations (Spearman) and identifies candidate time-lagged associations. Results go to `temporal_lagged_correlations.csv`.
5. Each causal method script (`causal_discovery_dbn.py`, `causal_discovery_granger.py`, `causal_discovery_tigramite.py`) performs method-specific discovery across taxa/time series, producing the corresponding CSV outputs in results.
6. `causal_consensus.py` ingests method-specific outputs and produces `causal_consensus_edges.csv`, a prioritized list of consensus-directed edges with supporting evidence counts and method-level metadata.

Outputs and interpretation
--------------------------
- Method outputs (CSV rows) typically include:
  - source_taxon, target_taxon, lag (if applicable), score/p-value, method
- `causal_consensus_edges.csv` summarizes edges supported by multiple methods. Columns include:
  - `source`, `target`, `supported_by` (comma-separated methods), `consensus_score` (aggregate ranking), `median_lag`, `notes`
- Use consensus edges as hypotheses for directed microbe→microbe or microbe→clinical-state interactions. These are not proofs of biological causation but prioritized leads for experimental validation.
- `differential_abundance_results.csv` helps interpret which taxa change with disease state and whether they appear in causal edges.
- Temporal lag files help identify whether changes in one taxon precede changes in another (useful for causal directionality).

Typical commands
----------------
Examples to run core steps (run from repository root):

- Activate environment (Windows PowerShell):
  - .venv\Scripts\Activate.ps1
- Preprocess:
  - python preprocess_pipeline.py --input hmp2_ibd_metagenomics_cleaned.csv --output hmp2_ibd_metagenomics_preprocessed.csv
- Run statistical analysis:
  - python statistical_analysis.py --input hmp2_ibd_metagenomics_preprocessed.csv --output differential_abundance_results.csv
- Run temporal analysis:
  - python temporal_analysis.py --input hmp2_ibd_metagenomics_temporal.csv --output temporal_lagged_correlations.csv
- Run causal methods (example; each script may have its own flags):
  - python causal_discovery_dbn.py --input hmp2_ibd_metagenomics_temporal.csv --output causal_dbn_results.csv
  - python causal_discovery_granger.py --input hmp2_ibd_metagenomics_temporal.csv --output causal_granger_results.csv
  - python causal_discovery_tigramite.py --input hmp2_ibd_metagenomics_temporal.csv --output causal_tigramite_results.csv
- Aggregate consensus:
  - python causal_consensus.py --inputs results/causal_dbn_results.csv,results/causal_granger_results.csv,results/causal_tigramite_results.csv --output causal_consensus_edges.csv

Configuration and reproducibility
--------------------------------
- If scripts accept configuration files or CLI flags, set seeds for deterministic steps (e.g., random forest-based ranking) and record package versions (e.g., `pip freeze > requirements.txt`).
- Store intermediate preprocessed data in data and final generated CSVs in results so analysis can be re-run from any stage.
- Consider adding a `Dockerfile` or `environment.yml` for fully reproducible environments.

Best practices & caveats
------------------------
- Large raw tables are included here; avoid committing ever-larger binary artifacts. GitHub's file size limit is 100 MB per file — current data files are under that, but be mindful when adding new datasets.
- Causal discovery methods have assumptions (stationarity, sampling frequency, adequate timepoints). Carefully validate and report those assumptions in the manuscript.
- Use multiple methods and consensus ranking to reduce method-specific artifacts, and annotate edges with method-specific caveats.

What to cite / references
-------------------------
- Cite original data sources (HMP2 / IBD multi-omics cohort) and the packages/algorithms used (Tigramite, Granger causality, DBN approaches, etc.) when publishing results.

Author & contact
----------------
- Primary: Aditya Akshat Singh (repo owner)
- For questions or reproducibility requests, open an issue or email the project owner

