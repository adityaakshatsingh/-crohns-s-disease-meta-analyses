"""
Causal Discovery: Granger Causality Analysis (VAR Models)
=========================================================

Author: Aditya Akshat Singh (CSE Student) & Antigravity (AI Coding Assistant)
Date: June 2026
Description: This script executes the first part of Phase 5: Causal Discovery.
             1. Fits bivariate Vector Autoregressive (VAR) models of lag 1 for
                candidate species and Fecal Calprotectin across 58 participants.
             2. Performs Granger Causality testing in both directions:
                - Forward: Does species CLR Granger-cause Calprotectin?
                - Backward: Does Calprotectin Granger-cause species CLR?
             3. Combines p-values across the cohort using Fisher's meta-analysis method.
             4. Applies Benjamini-Hochberg FDR correction.
             5. Outputs results to results/causal_granger_results.csv.

Dependencies: pandas, numpy, scipy, statsmodels
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import chi2
from statsmodels.tsa.stattools import grangercausalitytests
from statsmodels.stats.multitest import multipletests
import warnings

# File paths
PREPROCESSED_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_preprocessed.csv"
OUTPUT_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\results\causal_granger_results.csv"
RESULTS_DIR = r"C:\Users\voltt\OneDrive\Desktop\heal\results"

# Ignore the statsmodels verbose deprecation warning
warnings.filterwarnings('ignore', category=FutureWarning)

def run_granger_causality():
    print("------------------------------------------------------------")
    print("Starting Phase 5: Granger Causality Analysis (VAR)...")
    print("------------------------------------------------------------")

    if not os.path.exists(PREPROCESSED_FILE):
        print(f"Error: Preprocessed file not found at: {PREPROCESSED_FILE}")
        return

    # Load dataset
    df = pd.read_csv(PREPROCESSED_FILE)
    print(f"Loaded {df.shape[0]} samples.")

    # Select candidate species representing key findings from Phase 3 & 4
    candidates = [
        'Lawsonibacter asaccharolyticus_CLR',
        'Escherichia coli_CLR',
        'Haemophilus parainfluenzae_CLR',
        'Phascolarctobacterium faecium_CLR',
        'Dialister invisus_CLR',
        'Bifidobacterium longum_CLR',
        'Alistipes finegoldii_CLR',
        'Akkermansia muciniphila_CLR',
        'Paraprevotella xylaniphila_CLR',
        'Faecalibacterium prausnitzii_CLR',
        'Roseburia faecis_CLR'
    ]

    # Filter to only keep columns that are in the dataset
    candidates = [c for c in candidates if c in df.columns]
    print(f"Testing Granger Causality for {len(candidates)} candidate species.")

    # Group by participant
    grouped = df.groupby('Participant ID')
    
    results = []

    for species in candidates:
        p_vals_forward = []  # Species -> Calprotectin
        p_vals_backward = [] # Calprotectin -> Species
        
        for pid, group in grouped:
            group_sorted = group.sort_values(by='week_num')
            
            # Drop rows with missing calprotectin for VAR fitting
            group_clean = group_sorted.dropna(subset=['fecalcal'])
            
            # A VAR(1) model requires at least 6 points to have stable degrees of freedom
            if len(group_clean) < 6:
                continue
                
            data_var = group_clean[['fecalcal', species]].values
            
            # Check for zero variance
            if np.std(data_var[:, 0]) < 1e-6 or np.std(data_var[:, 1]) < 1e-6:
                continue
                
            # 1. Forward test: Species -> Calprotectin
            try:
                # statsmodels: 1st col is target, 2nd is predictor
                res_forward = grangercausalitytests(data_var, maxlag=1, verbose=False)
                # ssr_ftest returns (F-statistic, p-value, df_denom, df_num)
                p_f = res_forward[1][0]['ssr_ftest'][1]
                p_vals_forward.append(p_f)
            except:
                pass
                
            # 2. Backward test: Calprotectin -> Species
            data_var_rev = group_clean[[species, 'fecalcal']].values
            try:
                res_backward = grangercausalitytests(data_var_rev, maxlag=1, verbose=False)
                p_b = res_backward[1][0]['ssr_ftest'][1]
                p_vals_backward.append(p_b)
            except:
                pass
                
        # Combine p-values using Fisher's method
        def fisher_combine(p_list):
            p_list = [p for p in p_list if p is not None and not np.isnan(p)]
            k = len(p_list)
            if k == 0:
                return np.nan, np.nan
            # Clamp p-values to prevent log(0)
            p_list = [max(p, 1e-15) for p in p_list]
            stat = -2.0 * sum(np.log(p_list))
            combined_p = chi2.sf(stat, 2 * k)
            return stat, combined_p
            
        stat_f, p_f_combined = fisher_combine(p_vals_forward)
        stat_b, p_b_combined = fisher_combine(p_vals_backward)
        
        results.append({
            'Species': species.replace('_CLR', ''),
            'Forward_Fisher_Stat': stat_f,
            'Forward_p_value': p_f_combined,
            'Forward_N_Participants': len(p_vals_forward),
            'Backward_Fisher_Stat': stat_b,
            'Backward_p_value': p_b_combined,
            'Backward_N_Participants': len(p_vals_backward)
        })

    results_df = pd.DataFrame(results)

    # Apply FDR correction on forward and backward paths
    valid_f = results_df.dropna(subset=['Forward_p_value'])
    _, fdr_f, _, _ = multipletests(valid_f['Forward_p_value'], alpha=0.05, method='fdr_bh')
    results_df.loc[valid_f.index, 'Forward_FDR_p'] = fdr_f

    valid_b = results_df.dropna(subset=['Backward_p_value'])
    _, fdr_b, _, _ = multipletests(valid_b['Backward_p_value'], alpha=0.05, method='fdr_bh')
    results_df.loc[valid_b.index, 'Backward_FDR_p'] = fdr_b

    # Save to CSV
    os.makedirs(RESULTS_DIR, exist_ok=True)
    results_df.to_csv(OUTPUT_FILE, index=False)
    print(f"   - Saved Granger Causality results to: {OUTPUT_FILE}")

    # Print summary
    print("\n------------------------------------------------------------")
    print("Granger Causality Causal Pathways (FDR < 0.05):")
    print("------------------------------------------------------------")
    for _, row in results_df.iterrows():
        sp = row['Species']
        f_sig = row['Forward_FDR_p'] < 0.05
        b_sig = row['Backward_FDR_p'] < 0.05
        
        if f_sig and b_sig:
            print(f"   * {sp} <---> Fecal Calprotectin [Bi-directional Feedback Loop]")
            print(f"     (Forward p={row['Forward_FDR_p']:.2e}, Backward p={row['Backward_FDR_p']:.2e})")
        elif f_sig:
            print(f"   * {sp} ---> Fecal Calprotectin [Causal Driver]")
            print(f"     (Forward p={row['Forward_FDR_p']:.2e})")
        elif b_sig:
            print(f"   * Fecal Calprotectin ---> {sp} [Causal Responder]")
            print(f"     (Backward p={row['Backward_FDR_p']:.2e})")
    print("------------------------------------------------------------")

if __name__ == "__main__":
    run_granger_causality()
