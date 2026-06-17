# Phase 3: Static Statistical Analysis and Differential Abundance Report

**Project:** Temporal Causal Discovery of Microbial Signatures Preceding Inflammatory Flares in Crohn's Disease  
**Stage:** Phase 3 (Statistical Analysis)  
**Status:** Completed  

---

## 1. Objectives & Rationale

Causal discovery requires a solid understanding of the static association structure of the dataset. The objective of Phase 3 was to identify significant microbial associations with specific cohorts (Crohn's Disease, Ulcerative Colitis, Healthy controls) and intestinal inflammation states (continuous Calprotectin levels and active flare vs. remission states).

This phase addressed three critical areas of inquiry:
1. **Taxonomic Alterations in IBD:** Identifying species that differentiate disease states from healthy baselines.
2. **Inflammatory Covariation:** Isolating species that correlate with the gut inflammation marker, Fecal Calprotectin.
3. **Flare Biomarkers:** Finding species differentially abundant during active flares compared to remission.

---

## 2. Methodology & Implementation

The analysis was executed via the script `src/statistical_analysis.py` using the following statistical procedures:

1. **Cohort Contrast (Mann-Whitney U Test):**
   - Performed non-parametric two-sided Mann-Whitney U tests to compare CLR-abundance values across CD vs. Healthy, UC vs. Healthy, and CD vs. UC groups.
2. **Spearman Rank Correlation:**
   - Evaluated the monotonic correlation between microbial CLR-abundance and Fecal Calprotectin. Spearman's $\rho$ was chosen to capture non-linear relationship shapes.
3. **Flare vs. Remission Categorization:**
   - Grouped CD and UC patient samples into **Active Flare** ($\text{Calprotectin} \ge 150\ \mu\text{g/g}$, $n=370$ samples) and **Remission** ($< 150\ \mu\text{g/g}$, $n=566$ samples).
   - Ran two-sided Mann-Whitney U tests between these subgroups.
4. **Cliff's Delta Effect Size ($d$):**
   - Calculated Cliff's Delta to quantify the magnitude of differences between groups:
     $$d = \frac{2U}{n_1 n_2} - 1$$
     Effect size boundaries: $|d| < 0.147$ (negligible), $|d| < 0.33$ (small), $|d| < 0.474$ (medium), and $|d| \ge 0.474$ (large).
5. **False Discovery Rate (FDR) Correction:**
   - Multi-hypothesis testing correction was performed using the **Benjamini-Hochberg (BH)** procedure. Significance was defined at an adjusted $p < 0.05$.

---

## 3. Results & Outputs

Detailed results are saved in three outputs within the `results/` folder:
- **`results/differential_abundance_results.csv`**: Cohort comparisons.
- **`results/spearman_correlation_results.csv`**: Correlation coefficients.
- **`results/flare_vs_remission_results.csv`**: Subgroup analysis.

### Top Significant Findings (FDR-adjusted $p < 0.05$)

#### A. Crohn's Disease Dysbiosis (87 Significant Species)
- *Paraprevotella xylaniphila*: Enriched in CD (FDR $p=4.84 \times 10^{-27}$, $d=+0.484$, Large effect)
- *Clostridium sp. CAG 242*: Enriched in CD (FDR $p=4.84 \times 10^{-27}$, $d=+0.481$, Large effect)
- *Erysipelatoclostridium ramosum*: Enriched in CD (FDR $p=1.09 \times 10^{-24}$, $d=+0.459$, Medium effect)

#### B. Ulcerative Colitis Dysbiosis (68 Significant Species)
- *Alistipes finegoldii*: Depleted in UC (FDR $p=2.12 \times 10^{-18}$, $d=-0.453$, Medium effect)
- *Akkermansia muciniphila*: Depleted in UC (FDR $p=1.91 \times 10^{-11}$, $d=-0.351$, Medium effect)

#### C. Gut Inflammation (Calprotectin) Correlations (78 Significant Species)
- *Escherichia coli*: Positively correlated ($\rho = +0.272$, FDR $p=1.67 \times 10^{-19}$)
- *Lawsonibacter asaccharolyticus*: Negatively correlated ($\rho = -0.267$, FDR $p=3.26 \times 10^{-19}$)

#### D. Active Flares vs. Remission Biomarkers (55 Significant Species)
- *Lawsonibacter asaccharolyticus*: Depleted in flares (FDR $p=2.03 \times 10^{-15}$, $d=-0.329$)
- *Dialister invisus*: Enriched in flares (FDR $p=2.62 \times 10^{-14}$, $d=+0.312$)

---

## 4. Research-Grade Conclusions & Next-Step Preparedness

### Pathological Insights
- **The Enterobacteriaceae Signature:** The strong positive correlation of *Escherichia coli* with Calprotectin and its significant enrichment during active flares ($d = +0.277$) points to an inflammatory microenvironment. Gut inflammation releases oxygen and nitrate, promoting the expansion of facultative anaerobes (like *E. coli*) and displacing beneficial obligate anaerobes.
- **Mucosal Barrier Degradation in UC:** The significant depletion of *Akkermansia muciniphila* ($d = -0.351$) in Ulcerative Colitis represents a critical mucosal barrier signature. *A. muciniphila* degrades mucin to stimulate goblet cells and maintain the integrity of the epithelial lining; its depletion aligns with UC pathophysiology.
- **Inflammatory Depletion of SCFA Producers:** The depletion of *Lawsonibacter asaccharolyticus* during active flares ($d = -0.329$) suggests a loss of protective commensals, leading to reduced short-chain fatty acid (SCFA) production and exacerbation of gut inflammation.

These static associations provide a biological foundation. In the next phase, we transition from static correlation to temporal sequence mapping.
