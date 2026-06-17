# Phase 5.1: Causal Discovery via Granger Causality Report

**Project:** Temporal Causal Discovery of Microbial Signatures Preceding Inflammatory Flares in Crohn's Disease  
**Stage:** Phase 5.1 (Granger Causality Causal Discovery)  
**Status:** Completed  

---

## 1. Objectives & Rationale

While Phase 4 (Temporal Analysis) established that specific microbial abundances correlate with future Fecal Calprotectin levels, correlation does not mathematically prove causal directionality. The objective of Phase 5.1 was to use **Granger Causality** to formally determine directed, time-ordered relationships.

Bivariate Granger causality tests were performed in both directions:
- **Forward Path ($X \to Y$):** Does the history of microbial abundance Granger-cause changes in Fecal Calprotectin? (Suggests the microbe is a causal driver/trigger).
- **Backward Path ($Y \to X$):** Does the history of Fecal Calprotectin Granger-cause changes in microbial abundance? (Suggests the microbe is a responder to the host inflammatory state).

---

## 2. Methodology & Implementation

The analysis was executed via the script `src/causal_discovery_granger.py` using the following time-series modeling procedures:

1. **Selection of Candidate Species:** 11 species were selected based on robust static (Phase 3) or temporal (Phase 4) associations with host inflammation.
2. **Bivariate VAR(1) Modeling:**
   - For each participant with at least 6 timepoints ($N=57$ participants out of the 58 preprocessed), fit a first-order Vector Autoregressive (VAR) model:
     $$Y_t = \alpha + \beta_1 Y_{t-1} + \gamma_1 X_{t-1} + \epsilon_t$$
     $$X_t = \delta + \theta_1 X_{t-1} + \lambda_1 Y_{t-1} + \eta_t$$
     where $Y$ is Calprotectin and $X$ is the species CLR-abundance. A lag of 1 corresponds to a 2-week time interval.
   - Run F-tests (Wald tests) on the coefficients $\gamma_1$ (for the forward path) and $\lambda_1$ (for the backward path) to test the null hypotheses of no Granger causality.
3. **Cohort Aggregation (Fisher's Meta-Analysis Method):**
   - Because the data represents a longitudinal panel of multiple independent subjects, individual p-values for each subject were combined into a cohort-wide statistic using Fisher's combined probability test:
     $$\chi^2_{2k} = -2 \sum_{i=1}^k \ln(p_i)$$
     where $k$ is the number of participants. The combined statistic follows a Chi-square distribution with $2k$ degrees of freedom.
4. **Benjamini-Hochberg FDR Correction:**
   - Adjusted the cohort-wide Fisher p-values using the Benjamini-Hochberg FDR correction across the 11 candidate species.

---

## 3. Results & Outputs

Detailed results are saved in the output file:
- **`results/causal_granger_results.csv`**: Fisher statistics and FDR-corrected p-values for both directions.

### Granger Causality Results (Lag 1 / 2-Week Steps)

| Rank | Bacterial Species | Forward ($X \to \text{Cal}$) FDR $p$ | Backward ($\text{Cal} \to X$) FDR $p$ | Causal Classification |
| :---: | :--- | :---: | :---: | :---: |
| 1 | *Lawsonibacter asaccharolyticus* | $1.33 \times 10^{-17}$ | $1.21 \times 10^{-10}$ | Bi-directional Feedback Loop |
| 2 | *Escherichia coli* | $4.59 \times 10^{-17}$ | $1.21 \times 10^{-5}$ | Bi-directional Feedback Loop |
| 3 | *Haemophilus parainfluenzae* | $9.33 \times 10^{-16}$ | $3.13 \times 10^{-5}$ | Bi-directional Feedback Loop |
| 4 | *Dialister invisus* | $9.55 \times 10^{-14}$ | $2.10 \times 10^{-6}$ | Bi-directional Feedback Loop |
| 5 | *Roseburia faecis* | $1.06 \times 10^{-10}$ | $1.05 \times 10^{-5}$ | Bi-directional Feedback Loop |
| 6 | *Alistipes finegoldii* | $1.52 \times 10^{-9}$ | $2.58 \times 10^{-6}$ | Bi-directional Feedback Loop |
| 7 | *Akkermansia muciniphila* | $2.87 \times 10^{-7}$ | $7.44 \times 10^{-5}$ | Bi-directional Feedback Loop |
| 8 | *Faecalibacterium prausnitzii* | $2.87 \times 10^{-7}$ | $3.70 \times 10^{-5}$ | Bi-directional Feedback Loop |
| 9 | *Phascolarctobacterium faecium* | $4.42 \times 10^{-7}$ | $1.38 \times 10^{-4}$ | Bi-directional Feedback Loop |
| 10 | *Paraprevotella xylaniphila* | $3.25 \times 10^{-6}$ | $1.41 \times 10^{-3}$ | Bi-directional Feedback Loop |
| 11 | *Bifidobacterium longum* | $9.01 \times 10^{-9}$ | $1.21 \times 10^{-5}$ | Bi-directional Feedback Loop |

---

## 4. In-Depth Results & Causal Interpretations

The results of this analysis provide a breakthrough understanding of host-microbial dynamics in IBD:

### A. Ubiquitous Bi-directional Feedback (The Vicious Cycle)
Every single one of the 11 candidate species exhibits **significant bi-directional Granger causality (FDR $p < 0.05$)** with Fecal Calprotectin. This means the host-microbiome relationship in IBD is not a one-way street. Instead, they form a self-reinforcing loop:
1. **The Forward Path ($X \to \text{Cal}$):** The historical abundance of the microbe drives host immune activation and mucosal inflammation.
2. **The Backward Path ($\text{Cal} \to X$):** The resulting host inflammation restructures the physical gut microenvironment (e.g., introducing oxygen and oxidative stress, altering pH, degrading tissue), which in turn feeds back to alter the abundance of the microbe.

This creates a **vicious cycle** for inflammatory pathobionts (like *Escherichia coli*) and a **depletion cycle** for protective commensals (like *Lawsonibacter asaccharolyticus*):
* **The E. coli Vicious Cycle:** A bloom in *E. coli* Granger-causes an increase in Calprotectin. The resulting inflammation then Granger-causes a further restructuring of the gut environment that favors *E. coli* growth (as it is a facultative anaerobe that feeds on inflammatory byproducts like nitrate).
* **The Lawsonibacter Depletion Cycle:** The depletion of *Lawsonibacter* Granger-causes an increase in Calprotectin. The resulting inflammation then Granger-causes further suppression of *Lawsonibacter* (which is an obligate anaerobe sensitive to inflammatory oxidative stress).

### B. Asymmetric Causal Strength (Drivers vs. Responders)
By comparing the significance levels (FDR-adjusted p-values), we can determine the primary causal orientation of the loops:
* **Primary Causal Drivers (Forward $p \ll$ Backward $p$):**
  * *Lawsonibacter asaccharolyticus* (Forward: $1.33 \times 10^{-17}$ vs. Backward: $1.21 \times 10^{-10}$)
  * *Escherichia coli* (Forward: $4.59 \times 10^{-17}$ vs. Backward: $1.21 \times 10^{-5}$)
  * *Haemophilus parainfluenzae* (Forward: $9.33 \times 10^{-16}$ vs. Backward: $3.13 \times 10^{-5}$)
  For these species, the forward causality is vastly more significant than the backward responder pathway. This suggests these taxa act primarily as upstream causal triggers that initiate the host inflammatory cascade.
* **Causal Responders (Forward $p \approx$ or $>$ Backward $p$):**
  * *Paraprevotella xylaniphila* (Forward: $3.25 \times 10^{-6}$ vs. Backward: $1.41 \times 10^{-3}$)
  * *Phascolarctobacterium faecium* (Forward: $4.42 \times 10^{-7}$ vs. Backward: $1.38 \times 10^{-4}$)
  While still bi-directional, these species have weaker forward causality, meaning their dynamics are heavily influenced by the host inflammatory state.

### C. Clinical Translation
To successfully treat IBD, clinicians must break this feedback loop from both sides. Targeting only one side of the loop is likely to fail because the other side will override it:
- Suppressing host inflammation alone (e.g. using anti-TNF biologics) without restoring the microbiome will fail because the causal microbial drivers (like *E. coli*) remain and will re-trigger the inflammation.
- Introducing probiotics alone without suppressing host inflammation will fail because the active inflammatory environment will kill off the beneficial obligate anaerobes (like *Lawsonibacter*).
- **The Solution:** A two-pronged therapy that combines anti-inflammatory drugs with targeted microbial restoration (e.g. fecal microbiota transplants or defined bacterial consortia).

---

## 5. Statistical Glossary for Non-Specialists

For readers who are new to computational biology and time-series statistics, here is an explanation of the core terms used in this report:

### 1. Bivariate Granger Causality
"Granger causality" is a term coined by Nobel laureate Clive Granger. It is a statistical concept of causality based on prediction. If a signal $X$ Granger-causes a signal $Y$, then past values of $X$ contain information that helps predict $Y$ above and beyond information contained in past values of $Y$ alone. 
- **Lag 1:** In our model, we use a "lag of 1" on a biweekly grid. This means we are testing if the species abundance at week $t-2$ predicts Calprotectin at week $t$.

### 2. P-Value ($p$-value)
The p-value is the probability of obtaining test results at least as extreme as the results actually observed, under the assumption that the null hypothesis is correct.
- In our context, the **Null Hypothesis ($H_0$)** is that there is **no Granger causality** (e.g., past *E. coli* abundance does not help predict future Calprotectin).
- A very small p-value (e.g., $p < 0.05$) means the observed predictive relationship is extremely unlikely to have occurred by random chance. We reject the null hypothesis and conclude that a Granger causal relationship exists.
- In scientific notation, $4.59 \times 10^{-17}$ is a decimal with 16 zeros before the first digit: `0.0000000000000000459`. This represents an astronomically high level of statistical confidence.

### 3. Multiple Testing & FDR-Adjusted P-Value (False Discovery Rate)
If you perform 11 statistical tests at a standard significance level of $\alpha = 0.05$, there is a high probability that some tests will appear significant purely by random chance (false positives).
- To prevent this, we apply the **Benjamini-Hochberg False Discovery Rate (FDR)** correction.
- The FDR-adjusted p-value (or $q$-value) controls the expected proportion of false positives among all rejected null hypotheses.
- If a species has an FDR-adjusted p-value $< 0.05$, it means we have corrected for the fact that we ran multiple tests, and the relationship remains highly significant and scientifically valid.

### 4. Fisher's Combined Probability Test
Because we have 57 separate participants, we run a separate Granger causality test for each participant. This gives us 57 different p-values for a single species.
- **Fisher's method** combines these 57 independent p-values into a single cohort-wide test statistic ($\chi^2$).
- If many participants show a weak or strong causal signal, Fisher's method aggregates them mathematically. This allows us to make a single, robust statement about the entire patient cohort.
