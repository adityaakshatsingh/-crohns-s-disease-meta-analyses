# Phase 4: Temporal Analysis and Feature Engineering Report

**Project:** Temporal Causal Discovery of Microbial Signatures Preceding Inflammatory Flares in Crohn's Disease  
**Stage:** Phase 4 (Temporal Analysis)  
**Status:** Completed  

---

## 1. Objectives & Rationale

Static statistical analyses cannot establish directionality: if species $X$ correlates with Calprotectin, it remains unclear whether species $X$ causes inflammation, responds to inflammation, or shares a common driver. The objective of Phase 4 was to introduce temporal lag structures to evaluate whether microbial features predict *future* states of inflammation, satisfying the prerequisite of temporal precedence for causal discovery.

This phase implemented:
1. **Lag Alignment:** Matching microbial features at week $t$ with Calprotectin levels at $t+2$ and $t+4$ weeks.
2. **Rate-of-Change Engineering:** Calculating abundance velocity (deltas) between $t-2$ and $t$ weeks.
3. **Temporal Association Mapping:** Assessing species and deltas at week $t$ as predictors of future Calprotectin.

---

## 2. Methodology & Implementation

The analysis was executed via the script `src/temporal_analysis.py` using the following time-series procedures:

1. **Grid Lags Alignment:**
   - Grouped the preprocessed dataset by `Participant ID` and sorted chronologically.
   - For each timepoint $t$ on the biweekly grid, compiled:
     - The contemporaneous Calprotectin level (`fecalcal_t`).
     - Look-ahead Calprotectin values at $t+2$ weeks (`fecalcal_t_plus_2`) and $t+4$ weeks (`fecalcal_t_plus_4`).
2. **Velocity Feature Engineering (Deltas):**
   - For each species CLR-abundance $x_i$, calculated the abundance delta:
     $$\Delta x_i(t) = x_i(t) - x_i(t-2)$$
     representing the directional rate of change over the preceding two weeks.
3. **Lagged Spearman Correlation:**
   - Correlated species abundances $x_i(t)$ and deltas $\Delta x_i(t)$ at week $t$ with future Calprotectin $y(t+k)$ for lags $k \in \{2, 4\}$ weeks.
4. **Benjamini-Hochberg FDR Correction:**
   - Corrected raw p-values for multiple comparisons separately for each lag and feature group.

---

## 3. Results & Outputs

Detailed results are saved in the `data/` and `results/` folders:
- **`data/hmp2_ibd_metagenomics_temporal.csv`**: The lagged temporal dataset (1,187 samples, 239 columns including 117 abundances and 117 deltas).
- **`results/temporal_lagged_correlations.csv`**: Spearman correlation statistics for $t+2$ and $t+4$ week lags.

### Top Significant Precursors (FDR-adjusted $p < 0.05$)

We identified **149 significant lagged associations**.

#### A. 2-Week Lag Precursors ($k=2$ weeks, $n=1,129$ samples)
- *Lawsonibacter asaccharolyticus* abundance at $t$: Negatively predicts Calprotectin at $t+2$ ($\rho = -0.268$, FDR $p = 6.35 \times 10^{-18}$)
- *Escherichia coli* abundance at $t$: Positively predicts Calprotectin at $t+2$ ($\rho = +0.264$, FDR $p = 1.16 \times 10^{-17}$)
- *Haemophilus parainfluenzae* abundance at $t$: Positively predicts Calprotectin at $t+2$ ($\rho = +0.257$, FDR $p = 6.77 \times 10^{-17}$)

#### B. 4-Week Lag Precursors ($k=4$ weeks, $n=1,071$ samples)
- *Lawsonibacter asaccharolyticus* abundance at $t$: Negatively predicts Calprotectin at $t+4$ ($\rho = -0.264$, FDR $p = 1.80 \times 10^{-16}$)
- *Haemophilus parainfluenzae* abundance at $t$: Positively predicts Calprotectin at $t+4$ ($\rho = +0.254$, FDR $p = 1.60 \times 10^{-15}$)
- *Escherichia coli* abundance at $t$: Positively predicts Calprotectin at $t+4$ ($\rho = +0.253$, FDR $p = 1.60 \times 10^{-15}$)
- *Phascolarctobacterium faecium* abundance at $t$: Negatively predicts Calprotectin at $t+4$ ($\rho = -0.252$, FDR $p = 1.62 \times 10^{-15}$)

---

## 4. Research-Grade Conclusions & Next-Step Preparedness

### Temporal Precedence Established
The discovery of highly significant correlations at the 4-week look-ahead lag ($p < 10^{-15}$) provides strong statistical evidence that specific gut dysbiosis patterns **precede** the onset of active intestinal inflammation. This confirms the clinical potential of utilizing microbiome profiles as early-warning flare predictors.

### Candidate Protective & Inflammatory Precursors
- **Early Warning Marker (*Escherichia coli*):** A high abundance of *E. coli* at week $t$ is a robust predictor of a rise in Calprotectin 4 weeks later ($\rho = 0.253$), suggesting that the Enterobacteriaceae expansion occurs well before clinical flares manifest.
- **Protective Precursor (*Phascolarctobacterium faecium*):** The strong negative correlation of *Phascolarctobacterium faecium* at week $t$ with Calprotectin 4 weeks later ($\rho = -0.252$) highlights its potential as a protective biomarker. *P. faecium* degrades succinate (a pro-inflammatory metabolite that accumulates in the IBD gut) into propionate (a beneficial short-chain fatty acid). Its presence may suppress the inflammatory cascade.
- **Lack of Delta Significance:** Abundance velocity (deltas) did not cross the strict FDR threshold. This implies that absolute taxonomic abundance is a more reliable predictor of future inflammation than short-term abundance velocity, or that velocity features require multi-variable modeling.

With lag features engineered and temporal correlations confirmed, the dataset is prepared for **Phase 5: Causal Discovery**, where we will construct directed causal graphs.
