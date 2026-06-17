"""
Causal Discovery: PCMCI Analysis (Tigramite Package)
=====================================================

Author: Aditya Akshat Singh (CSE Student) & Antigravity (AI Coding Assistant)
Date: June 2026
Description: This script executes the third part of Phase 5: Causal Discovery.
             1. Compiles longitudinal panel data into a single array with nan-spacers.
             2. Initializes Tigramite DataFrame, mask, and ParCorr independence test.
             3. Runs PCMCI (PC1 parent selection + MCI conditional independence testing)
                for lags tau=1 (2 weeks) and tau=2 (4 weeks).
             4. Applies Benjamini-Hochberg FDR correction across all lagged links.
             5. Outputs results to results/causal_tigramite_results.csv.

Dependencies: pandas, numpy, scipy, statsmodels, tigramite
"""

import os
import numpy as np
import pandas as pd
import math
from statsmodels.stats.multitest import multipletests
import warnings

# Tigramite imports
import tigramite.data_processing as pp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.parcorr import ParCorr

# File paths
PREPROCESSED_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_preprocessed.csv"
OUTPUT_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\results\causal_tigramite_results.csv"
RESULTS_DIR = r"C:\Users\voltt\OneDrive\Desktop\heal\results"

# Ignore warnings
warnings.filterwarnings('ignore')

def run_pcmci_analysis():
    print("------------------------------------------------------------")
    print("Starting Phase 5.3: PCMCI Causal Discovery (Tigramite)...")
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
    print(f"Variables included in PCMCI: {len(variables)} nodes.")

    # Group by participant
    grouped = df.groupby('Participant ID')
    compiled_data = []

    for pid, group in grouped:
        group_sorted = group.sort_values(by='week_num')
        
        times = [float(row['week_num']) for _, row in group_sorted.iterrows()]
        min_week = min(times)
        max_week = max(times)
        
        # Determine biweekly grid
        start_grid = int(math.ceil(min_week / 2.0) * 2)
        end_grid = int(math.floor(max_week / 2.0) * 2)
        grid = list(range(start_grid, end_grid + 1, 2))
        
        # Map week_num -> row
        week_map = {float(row['week_num']): row for _, row in group_sorted.iterrows()}
        
        p_data = []
        for w in grid:
            w_float = float(w)
            row = week_map.get(w_float, None)
            
            row_vals = []
            for var in variables:
                if row is not None and not pd.isna(row[var]):
                    row_vals.append(float(row[var]))
                else:
                    row_vals.append(np.nan)
            p_data.append(row_vals)
            
        compiled_data.extend(p_data)
        # Append a separator row of np.nan to prevent cross-participant causal links
        compiled_data.append([np.nan] * len(variables))

    # Convert to numpy array (drop the final nan-spacer row)
    data_array = np.array(compiled_data[:-1])
    print(f"Compiled panel time-series shape: {data_array.shape}")

    var_names = [v.replace('_CLR', '') for v in variables]
    mask = np.isnan(data_array)

    # Fill NaNs with 0.0 in the data array, since mask tells Tigramite to ignore them
    data_array_filled = np.nan_to_num(data_array, nan=0.0)

    # Initialize Tigramite DataFrame
    dataframe = pp.DataFrame(data_array_filled, var_names=var_names, mask=mask)

    # Initialize ParCorr test (Partial Correlation) using analytic significance
    cond_ind_test = ParCorr(significance='analytic')

    # Initialize PCMCI
    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=cond_ind_test, verbosity=0)

    # Run PCMCI (PC1 + MCI steps) with tau_max = 2 (lags 1 and 2, corresponding to 2 and 4 weeks)
    print("Running PCMCI algorithm...")
    results = pcmci.run_pcmci(tau_max=2, pc_alpha=0.05, alpha_level=0.05)

    p_matrix = results['p_matrix']
    val_matrix = results['val_matrix']

    links = []
    num_nodes = len(var_names)

    for i in range(num_nodes):
        for j in range(num_nodes):
            for tau in [1, 2]:
                p_val = p_matrix[i, j, tau]
                val = val_matrix[i, j, tau]
                
                if not np.isnan(p_val):
                    links.append({
                        'Predictor': var_names[i],
                        'Target': var_names[j],
                        'Lag_Weeks': int(tau * 2),
                        'MCI_Correlation': val,
                        'p_value': p_val
                    })

    links_df = pd.DataFrame(links)

    # Apply FDR correction
    _, fdr_p, _, _ = multipletests(links_df['p_value'], alpha=0.05, method='fdr_bh')
    links_df['FDR_adjusted_p'] = fdr_p

    # Sort results
    links_df = links_df.sort_values(by='FDR_adjusted_p')

    # Save to CSV
    os.makedirs(RESULTS_DIR, exist_ok=True)
    links_df.to_csv(OUTPUT_FILE, index=False)
    print(f"   - Saved PCMCI results to: {OUTPUT_FILE}")

    # Print significant cross-variable edges
    print("\n------------------------------------------------------------")
    print("Significant Cross-Variable PCMCI Causal Edges (FDR < 0.05):")
    print("------------------------------------------------------------")
    cross_sig = links_df[(links_df['FDR_adjusted_p'] < 0.05) & (links_df['Predictor'] != links_df['Target'])]
    
    for _, row in cross_sig.iterrows():
        p_val = row['FDR_adjusted_p']
        mci_val = row['MCI_Correlation']
        pred = row['Predictor']
        tar = row['Target']
        lag = row['Lag_Weeks']
        
        rel_str = "positively drives (promotes)" if mci_val > 0 else "negatively drives (inhibits)"
        print(f"   * {pred} (t-{lag}) {rel_str} {tar} (t) [MCI={mci_val:.3f}, FDR p={p_val:.2e}]")
    print("------------------------------------------------------------")

if __name__ == "__main__":
    run_pcmci_analysis()
