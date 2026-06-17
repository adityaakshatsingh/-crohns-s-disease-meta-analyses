"""
Causal Discovery: Dynamic Bayesian Network (DBN) Analysis
=========================================================

Author: Aditya Akshat Singh (CSE Student) & Antigravity (AI Coding Assistant)
Date: June 2026
Description: This script executes the second part of Phase 5: Causal Discovery.
             1. Compiles pooled transition pairs (t-2 -> t) across all 58 participants.
             2. Fits multivariate Ordinary Least Squares (OLS) models to learn
                the transition parameters of the Dynamic Bayesian Network.
             3. Runs OLS regressions for each target node using all lagged variables
                simultaneously as features.
             4. Applies Benjamini-Hochberg FDR correction across all 144 edges.
             5. Outputs results to results/causal_dbn_results.csv.

Dependencies: pandas, numpy, scipy, statsmodels
"""

import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import warnings

# File paths
PREPROCESSED_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_preprocessed.csv"
OUTPUT_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\results\causal_dbn_results.csv"
RESULTS_DIR = r"C:\Users\voltt\OneDrive\Desktop\heal\results"

# Ignore statsmodels warnings
warnings.filterwarnings('ignore')

def run_dbn_analysis():
    print("------------------------------------------------------------")
    print("Starting Phase 5.2: Dynamic Bayesian Network (DBN) Analysis...")
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
    variables = ['fecalcal'] + candidates
    print(f"Variables included in DBN: {len(variables)} nodes.")

    # Group by participant and construct pooled transition pairs (t-2 -> t)
    grouped = df.groupby('Participant ID')
    transition_data = []

    for pid, group in grouped:
        group_sorted = group.sort_values(by='week_num')
        
        # Drop rows where calprotectin is missing to keep clinical measurements clean
        group_clean = group_sorted.dropna(subset=['fecalcal'])
        
        # Map week_num -> row
        week_map = {float(row['week_num']): row for _, row in group_clean.iterrows()}
        weeks = sorted(list(week_map.keys()))
        
        for t in weeks:
            if (t - 2.0) in week_map:
                row_t_minus_2 = week_map[t - 2.0]
                row_t = week_map[t]
                
                # Predictors at t-2 (Lag 1)
                predictors = {}
                for var in variables:
                    predictors[f"{var}_lag1"] = float(row_t_minus_2[var])
                    
                # Targets at t (Current)
                targets = {}
                for var in variables:
                    targets[f"{var}_t"] = float(row_t[var])
                    
                # Combine
                pair = {**predictors, **targets}
                transition_data.append(pair)

    transition_df = pd.DataFrame(transition_data)
    print(f"Compiled {transition_df.shape[0]} pooled time transition pairs.")

    # Run OLS for each target variable using all lagged variables as features
    dbn_results = []

    for target_var in variables:
        y_col = f"{target_var}_t"
        X_cols = [f"{var}_lag1" for var in variables]
        
        y = transition_df[y_col].values
        X = transition_df[X_cols].values
        
        # Add intercept to regression
        X_const = sm.add_constant(X)
        
        # Fit OLS
        model = sm.OLS(y, X_const)
        model_fit = model.fit()
        
        # Extract parameters, t-statistics, and p-values
        # Index 0 is constant, 1 onwards are predictors
        params = model_fit.params[1:]
        tvalues = model_fit.tvalues[1:]
        pvalues = model_fit.pvalues[1:]
        
        for idx, pred_var in enumerate(variables):
            dbn_results.append({
                'Target': target_var.replace('_CLR', ''),
                'Predictor': pred_var.replace('_CLR', ''),
                'Coefficient': params[idx],
                't_Statistic': tvalues[idx],
                'p_value': pvalues[idx]
            })

    dbn_df = pd.DataFrame(dbn_results)

    # Apply FDR correction across all 144 equations (12 targets * 12 predictors)
    _, fdr_p, _, _ = multipletests(dbn_df['p_value'], alpha=0.05, method='fdr_bh')
    dbn_df['FDR_adjusted_p'] = fdr_p

    # Save to CSV
    os.makedirs(RESULTS_DIR, exist_ok=True)
    dbn_df.to_csv(OUTPUT_FILE, index=False)
    print(f"   - Saved DBN results to: {OUTPUT_FILE}")

    # Print significant cross-variable edges (excluding self-loops for clarity)
    print("\n------------------------------------------------------------")
    print("Significant Cross-Variable DBN Edges (FDR < 0.05):")
    print("------------------------------------------------------------")
    cross_sig = dbn_df[(dbn_df['FDR_adjusted_p'] < 0.05) & (dbn_df['Predictor'] != dbn_df['Target'])]
    cross_sig_sorted = cross_sig.sort_values(by='FDR_adjusted_p')
    
    for _, row in cross_sig_sorted.iterrows():
        p_val = row['FDR_adjusted_p']
        coef = row['Coefficient']
        pred = row['Predictor']
        tar = row['Target']
        
        rel_str = "promotes" if coef > 0 else "inhibits"
        print(f"   * {pred} (t-2) {rel_str} {tar} (t) [Coef={coef:.3f}, FDR p={p_val:.2e}]")
    print("------------------------------------------------------------")

if __name__ == "__main__":
    run_dbn_analysis()
