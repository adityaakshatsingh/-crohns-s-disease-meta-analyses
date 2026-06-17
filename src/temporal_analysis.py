"""
Temporal Analysis and Feature Engineering for IBD Metagenomics Dataset
======================================================================

Author: Aditya Akshat Singh (CSE Student) & Antigravity (AI Coding Assistant)
Date: June 2026
Description: This script executes Phase 4 of the research plan:
             1. Resamples and aligns time-series data to build 2-week and 4-week lags.
             2. Engineers abundance deltas (change in relative abundance from t-2 to t).
             3. Runs Spearman correlations between features at week t and future 
                inflammation (Calprotectin) at weeks t+2 and t+4 (FDR-corrected).
             4. Generates a final model-ready temporal dataset.
             
Dependencies: pandas, numpy, scipy, statsmodels
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# File paths
PREPROCESSED_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_preprocessed.csv"
OUTPUT_TEMPORAL_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_temporal.csv"
RESULTS_DIR = r"C:\Users\voltt\OneDrive\Desktop\heal\results"

def run_temporal_analysis():
    print("------------------------------------------------------------")
    print("Starting Phase 4: Temporal Analysis & Feature Engineering...")
    print("------------------------------------------------------------")

    if not os.path.exists(PREPROCESSED_FILE):
        print(f"Error: Preprocessed file not found at: {PREPROCESSED_FILE}")
        return

    # Load dataset
    df = pd.read_csv(PREPROCESSED_FILE)
    print(f"Loaded {df.shape[0]} samples and {df.shape[1]} columns.")

    species_cols = [col for col in df.columns if col.endswith('_CLR')]

    # Group by participant to construct lags and deltas per individual
    grouped = df.groupby('Participant ID')
    
    temporal_rows = []

    for pid, group in grouped:
        # Sort chronologically by week_num
        group_sorted = group.sort_values(by='week_num').copy()
        
        # Build dictionary mapping week_num -> row for quick lookup
        week_map = {}
        for _, row in group_sorted.iterrows():
            week_map[float(row['week_num'])] = row
            
        weeks = sorted(list(week_map.keys()))
        
        for idx, t in enumerate(weeks):
            current_row = week_map[t]
            diagnosis = current_row['diagnosis']
            ext_id = current_row['External ID']
            
            fcal_t = float(current_row['fecalcal']) if not pd.isna(current_row['fecalcal']) else np.nan
            
            # 1. Look ahead for future Calprotectin values
            fcal_t2 = np.nan
            fcal_t4 = np.nan
            
            if (t + 2.0) in week_map:
                row_t2 = week_map[t + 2.0]
                fcal_t2 = float(row_t2['fecalcal']) if not pd.isna(row_t2['fecalcal']) else np.nan
                
            if (t + 4.0) in week_map:
                row_t4 = week_map[t + 4.0]
                fcal_t4 = float(row_t4['fecalcal']) if not pd.isna(row_t4['fecalcal']) else np.nan
                
            # 2. Look back to calculate Deltas (t minus t-2)
            has_prev = (t - 2.0) in week_map
            prev_row = week_map[t - 2.0] if has_prev else None
            
            deltas = {}
            for col in species_cols:
                current_val = float(current_row[col])
                if has_prev:
                    prev_val = float(prev_row[col])
                    deltas[f"{col}_delta"] = current_val - prev_val
                else:
                    deltas[f"{col}_delta"] = 0.0 # Default delta is 0 for boundary
                    
            # 3. Compile temporal features
            out_row = {
                'External ID': ext_id,
                'Participant ID': pid,
                'week_num': t,
                'diagnosis': diagnosis,
                'fecalcal_t': fcal_t,
                'fecalcal_t_plus_2': fcal_t2,
                'fecalcal_t_plus_4': fcal_t4
            }
            
            # Add species CLR abundances at t
            for col in species_cols:
                out_row[col] = float(current_row[col])
                
            # Add species CLR deltas at t
            for delta_col, val in deltas.items():
                out_row[delta_col] = val
                
            temporal_rows.append(out_row)

    temporal_df = pd.DataFrame(temporal_rows)
    temporal_df.to_csv(OUTPUT_TEMPORAL_FILE, index=False)
    print(f"   - Successfully saved temporal dataset to: {OUTPUT_TEMPORAL_FILE}")

    # ==========================================
    # Lagged Correlation Analysis
    # ==========================================
    print("\nRunning Lagged Correlation Analysis...")
    
    # We will correlate abundance at t with calprotectin at t+2 and t+4
    # Filter subsets where Calprotectin is available
    df_t2 = temporal_df.dropna(subset=['fecalcal_t_plus_2'])
    df_t4 = temporal_df.dropna(subset=['fecalcal_t_plus_4'])
    
    print(f"   - Samples available for t+2 lag analysis: {df_t2.shape[0]}")
    print(f"   - Samples available for t+4 lag analysis: {df_t4.shape[0]}")
    
    lagged_results = []
    
    # Perform Spearman correlations for t+2 lag
    if df_t2.shape[0] >= 5:
        fcal_vals = df_t2['fecalcal_t_plus_2'].values
        for col in species_cols:
            # 1. Abundance at t vs Calprotectin at t+2
            rho, p_val = spearmanr(df_t2[col].values, fcal_vals)
            # 2. Delta abundance at t vs Calprotectin at t+2
            rho_d, p_val_d = spearmanr(df_t2[f"{col}_delta"].values, fcal_vals)
            
            lagged_results.append({
                'Species': col.replace('_CLR', ''),
                'Lag': 't_plus_2',
                'Feature_Type': 'Abundance',
                'Spearman_Rho': rho,
                'p_value': p_val
            })
            lagged_results.append({
                'Species': col.replace('_CLR', ''),
                'Lag': 't_plus_2',
                'Feature_Type': 'Delta',
                'Spearman_Rho': rho_d,
                'p_value': p_val_d
            })
            
    # Perform Spearman correlations for t+4 lag
    if df_t4.shape[0] >= 5:
        fcal_vals = df_t4['fecalcal_t_plus_4'].values
        for col in species_cols:
            # 1. Abundance at t vs Calprotectin at t+4
            rho, p_val = spearmanr(df_t4[col].values, fcal_vals)
            # 2. Delta abundance at t vs Calprotectin at t+4
            rho_d, p_val_d = spearmanr(df_t4[f"{col}_delta"].values, fcal_vals)
            
            lagged_results.append({
                'Species': col.replace('_CLR', ''),
                'Lag': 't_plus_4',
                'Feature_Type': 'Abundance',
                'Spearman_Rho': rho,
                'p_value': p_val
            })
            lagged_results.append({
                'Species': col.replace('_CLR', ''),
                'Lag': 't_plus_4',
                'Feature_Type': 'Delta',
                'Spearman_Rho': rho_d,
                'p_value': p_val_d
            })

    if lagged_results:
        lagged_df = pd.DataFrame(lagged_results)
        
        # Apply FDR correction separately per lag & feature combination
        groups = lagged_df.groupby(['Lag', 'Feature_Type'])
        corrected_dfs = []
        for name, group in groups:
            grp_df = group.copy()
            _, fdr_p, _, _ = multipletests(grp_df['p_value'], alpha=0.05, method='fdr_bh')
            grp_df['FDR_adjusted_p'] = fdr_p
            corrected_dfs.append(grp_df)
            
        final_lagged_df = pd.concat(corrected_dfs, ignore_index=True)
        # Sort by significance
        final_lagged_df['Abs_Rho'] = final_lagged_df['Spearman_Rho'].abs()
        final_lagged_df = final_lagged_df.sort_values(by=['FDR_adjusted_p', 'Abs_Rho'], ascending=[True, False]).drop(columns=['Abs_Rho'])
        
        corr_output = os.path.join(RESULTS_DIR, 'temporal_lagged_correlations.csv')
        final_lagged_df.to_csv(corr_output, index=False)
        print(f"   - Saved lagged correlation analysis to: {corr_output}")
        
        # ==========================================
        # Report Discoveries
        # ==========================================
        print("\n------------------------------------------------------------")
        print("Top Temporal Precursor Discoveries (FDR < 0.05):")
        print("------------------------------------------------------------")
        
        sig_lagged = final_lagged_df[final_lagged_df['FDR_adjusted_p'] < 0.05]
        print(f"Found {len(sig_lagged)} significant lagged associations.")
        
        if not sig_lagged.empty:
            # Show top 8 significant lagged features
            for _, row in sig_lagged.head(10).iterrows():
                lag_str = "2 weeks later (t+2)" if row['Lag'] == 't_plus_2' else "4 weeks later (t+4)"
                feat_str = "absolute abundance" if row['Feature_Type'] == 'Abundance' else "rate of change (delta)"
                dir_str = "positively predicts" if row['Spearman_Rho'] > 0 else "negatively predicts (protective)"
                
                print(f"   * {row['Species']} {feat_str} {dir_str} inflammation {lag_str} (Rho={row['Spearman_Rho']:.3f}, FDR p={row['FDR_adjusted_p']:.4e})")
        else:
            print("No individual species abundances crossed the strict multiple-testing FDR < 0.05 significance threshold.")
            print("This suggests that predictive signals may require multi-variable models (such as machine learning classifiers) rather than univariate correlations.")
        print("------------------------------------------------------------")
    else:
        print("   - Error: Could not compute lagged correlation (insufficient samples).")

if __name__ == "__main__":
    run_temporal_analysis()
