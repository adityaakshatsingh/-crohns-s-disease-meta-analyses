# Phase 2: Data Preprocessing and Exploration Report

**Project:** Temporal Causal Discovery of Microbial Signatures Preceding Inflammatory Flares in Crohn's Disease  
**Stage:** Phase 2 (Preprocessing and Exploration)  
**Status:** Completed  

---

## 1. Objectives & Rationale

Metagenomic data derived from high-throughput sequencing presents unique statistical challenges, notably compositional constraints, high sparsity, longitudinal irregularity, and duplicate sequencing records. The objective of Phase 2 was to construct a robust, clean, and model-ready dataset that satisfies the assumptions of temporal causal discovery algorithms.

This phase resolved three major data quality issues:
1. **Redundancy & Duplication:** Deduplicating exact-matching rows and resolving conflicting records.
2. **Compositionality:** Transforming relative abundances to break mathematical constraints.
3. **Temporal Inconsistency:** Resampling irregular observation windows into a standardized biweekly grid.

---

## 2. Methodology & Implementation

The preprocessing pipeline was executed via the scripts `clean_dataset.py` and `preprocess_pipeline.py` using the following research-grade steps:

1. **Deduplication and Sample Merging:** 
   - Grouped the raw dataset by `External ID` (representing unique physical stool samples) to identify duplicate sequencing runs.
   - Retained the row with completed Fecal Calprotectin (`fecalcal`) data in cases of duplicate samples.
   - Sorted all samples chronologically by `Participant ID` and `week_num`.
   - Merged duplicate records for the same participant in the same week by averaging their microbial abundances and calprotectin values.
2. **Prevalence and Abundance Filtering:**
   - Filtered out rare taxa to reduce dimensionality and mitigate zero-inflation.
   - Retained species present in $>10\%$ of all samples with a mean relative abundance of $>0.01\%$.
3. **Participant Inclusion Criteria:**
   - Retained only participants with sufficient longitudinal data density.
   - Criteria: $\ge 4$ longitudinal stool samples and $\le 70\%$ missingness in Fecal Calprotectin.
4. **Biweekly Resampling and Linear Interpolation:**
   - For each kept participant, defined a regular target grid of even-numbered weeks (e.g., $0.0, 2.0, 4.0, \dots$).
   - Performed linear interpolation for both calprotectin and microbial relative abundances to align samples to this grid, using nearest-neighbor extrapolation at boundaries.
5. **Centered Log-Ratio (CLR) Transformation:**
   - Metagenomic relative abundances sum to 100%, causing a "closed composition" problem where variables are not independent.
   - Added a tiny pseudocount ($10^{-6}\%$) to replace zeros.
   - Applied the Centered Log-Ratio (CLR) transform to project the abundance data from the Simplex space to Euclidean space:
     $$\text{CLR}(x) = \left[ \ln\left(\frac{x_1}{g(x)}\right), \dots, \ln\left(\frac{x_D}{g(x)}\right) \right]$$
     where $g(x)$ is the geometric mean of the composition vector.

---

## 3. Results & Outputs

- **Raw Data Input:** 3,387 rows, 571 columns, 116 participants, 566 microbial species.
- **Cleaned Dataset:** 1,367 unique samples.
- **Filtered Feature Space:** Reduced microbial species columns from **566 to 117**.
- **Filtered Subject Cohort:** Retained **58 participants** (29 Crohn's Disease, 16 Ulcerative Colitis, 13 Healthy Controls) out of the original 116.
- **Final Model-Ready Dataset:** Generated the file `data/hmp2_ibd_metagenomics_preprocessed.csv` containing **1,187 samples** and **122 columns** (5 metadata columns + 117 CLR-transformed species abundances).

---

## 4. Research-Grade Conclusions & Next-Step Preparedness

### Compositional Rigor
By mapping compositional relative abundances into Euclidean space using the CLR transformation, we have mathematically eliminated the risk of *spurious correlation* (negative correlation bias). This is a vital prerequisite for causal discovery algorithms (like Tigramite/PCMCI or Dynamic Bayesian Networks), which assume baseline variable independence.

### Longitudinal Alignment
Resampling onto a biweekly grid resolves the irregular time steps of the raw dataset. This ensures that any time-lagged features engineered in future steps (e.g., $t-2$, $t-4$ weeks) correspond to fixed biological time steps, validating the assumptions of Granger Causality.

### Clinical Cohort Integrity
Filtering out participants with $>70\%$ missing calprotectin values ensures that our target variable (intestinal inflammation) is grounded in actual clinical measurements, rather than excessive interpolation. Retaining 58 highly dense participants (averaging ~20.5 weeks of longitudinal follow-up) provides a statistically robust sample size for training causal models.
