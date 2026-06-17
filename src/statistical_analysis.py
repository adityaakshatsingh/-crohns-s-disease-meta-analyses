"""
Research-Grade Statistical Analysis for HMP2 IBD Metagenomics Dataset
=====================================================================

Author: Aditya Akshat Singh (CSE Student) & Antigravity (AI Coding Assistant)
Date: June 2026
Description: This script performs differential abundance testing (Mann-Whitney U
             with Cliff's Delta effect size), Spearman rank correlation with Fecal 
             Calprotectin, and Active Flare vs. Remission differential analysis. 
             All p-values are adjusted using Benjamini-Hochberg FDR correction.
             
Dependencies: pandas, numpy, scipy, statsmodels
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

# File paths
PREPROCESSED_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_preprocessed.csv"
RESULTS_DIR = r"C:\Users\voltt\OneDrive\Desktop\heal\results"

def calculate_cliffs_delta(x, y):
    """
    Computes Cliff's Delta non-parametric effect size.
    d = (2 * U) / (n1 * n2) - 1
    Interpretation:
       |d| < 0.147 : Negligible
       |d| < 0.330 : Small
       |d| < 0.474 : Medium
       |d| >= 0.474: Large
    """
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return 0.0
    # mannwhitneyu returns U for the first sample (x)
    u_stat, _ = mannwhitneyu(x, y, alternative='two-sided')
    delta = (2.0 * u_stat) / (n1 * n2) - 1.0
    return delta

def run_analysis():
    print("------------------------------------------------------------")
    print("Starting Research-Grade Statistical Analysis Pipeline...")
    print("------------------------------------------------------------")

    if not os.path.exists(PREPROCESSED_FILE):
        print(f"Error: Preprocessed file not found at: {PREPROCESSED_FILE}")
        return

    # Create results directory if it doesn't exist
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Load dataset
    df = pd.read_csv(PREPROCESSED_FILE)
    print(f"Loaded {df.shape[0]} samples and {df.shape[1]} columns.")

    # Species columns are all columns ending with '_CLR'
    species_cols = [col for col in df.columns if col.endswith('_CLR')]
    print(f"Identified {len(species_cols)} preprocessed microbial species.")

    # Group cohorts
    cd_data = df[df['diagnosis'] == 'Crohns Disease']
    uc_data = df[df['diagnosis'] == 'Ulcerative Colitis']
    healthy_data = df[df['diagnosis'] == 'Healthy']

    print(f"Cohort sizes: CD = {cd_data.shape[0]} samples, UC = {uc_data.shape[0]} samples, Healthy = {healthy_data.shape[0]} samples.")

    # ==========================================
    # 1. Differential Abundance Analysis
    # ==========================================
    print("\n1. Running Differential Abundance Analysis (Mann-Whitney U)...")
    
    comparisons = [
        ('CD_vs_Healthy', cd_data, healthy_data),
        ('UC_vs_Healthy', uc_data, healthy_data),
        ('CD_vs_UC', cd_data, uc_data)
    ]

    diff_ab_results = []

    for comp_name, df1, df2 in comparisons:
        comp_results = []
        for col in species_cols:
            x = df1[col].values
            y = df2[col].values
            
            # Mann-Whitney U test
            u_stat, p_val = mannwhitneyu(x, y, alternative='two-sided')
            
            # Cliff's Delta effect size
            delta = calculate_cliffs_delta(x, y)
            
            # Means for reference
            mean1 = np.mean(x)
            mean2 = np.mean(y)
            mean_diff = mean1 - mean2
            
            comp_results.append({
                'Comparison': comp_name,
                'Species': col.replace('_CLR', ''),
                'Mean_Group1': mean1,
                'Mean_Group2': mean2,
                'Mean_Difference': mean_diff,
                'U_Statistic': u_stat,
                'p_value': p_val,
                'Cliffs_Delta': delta
            })
            
        # Apply Benjamini-Hochberg FDR correction for this comparison group
        comp_df = pd.DataFrame(comp_results)
        _, fdr_p, _, _ = multipletests(comp_df['p_value'], alpha=0.05, method='fdr_bh')
        comp_df['FDR_adjusted_p'] = fdr_p
        
        diff_ab_results.append(comp_df)

    final_diff_df = pd.concat(diff_ab_results, ignore_index=True)
    # Reorder columns
    cols_order = ['Comparison', 'Species', 'Mean_Group1', 'Mean_Group2', 'Mean_Difference', 'U_Statistic', 'p_value', 'FDR_adjusted_p', 'Cliffs_Delta']
    final_diff_df = final_diff_df[cols_order]
    
    diff_output = os.path.join(RESULTS_DIR, 'differential_abundance_results.csv')
    final_diff_df.to_csv(diff_output, index=False)
    print(f"   - Saved differential abundance results to: {diff_output}")

    # ==========================================
    # 2. Spearman Correlation with Calprotectin
    # ==========================================
    print("\n2. Running Spearman Rank Correlation with Fecal Calprotectin...")
    
    # Filter samples where Calprotectin is available
    df_fcal = df.dropna(subset=['fecalcal'])
    print(f"   - Found {df_fcal.shape[0]} samples with non-empty Calprotectin values.")
    
    corr_results = []
    fcal_vals = df_fcal['fecalcal'].values
    
    for col in species_cols:
        species_vals = df_fcal[col].values
        
        # Spearman correlation
        rho, p_val = spearmanr(species_vals, fcal_vals)
        
        corr_results.append({
            'Species': col.replace('_CLR', ''),
            'Spearman_Rho': rho,
            'p_value': p_val
        })
        
    corr_df = pd.DataFrame(corr_results)
    # Apply FDR correction
    _, fdr_p, _, _ = multipletests(corr_df['p_value'], alpha=0.05, method='fdr_bh')
    corr_df['FDR_adjusted_p'] = fdr_p
    
    # Sort by absolute correlation coefficient
    corr_df['Abs_Rho'] = corr_df['Spearman_Rho'].abs()
    corr_df = corr_df.sort_values(by='Abs_Rho', ascending=False).drop(columns=['Abs_Rho'])
    
    corr_output = os.path.join(RESULTS_DIR, 'spearman_correlation_results.csv')
    corr_df.to_csv(corr_output, index=False)
    print(f"   - Saved Spearman correlation results to: {corr_output}")

    # ==========================================
    # 3. Active Flare vs. Remission Subgroup Analysis
    # ==========================================
    print("\n3. Running Subgroup Analysis: Active Flare vs. Remission...")
    
    # Keep only patients with Crohn's or UC (exclude healthy controls for flare analysis)
    df_ibd = df_fcal[df_fcal['diagnosis'].isin(['Crohns Disease', 'Ulcerative Colitis'])].copy()
    
    # Define active flare vs. remission (clinical standard threshold: Calprotectin >= 150 ug/g)
    df_ibd['Flare_Status'] = df_ibd['fecalcal'].apply(lambda x: 'Active_Flare' if x >= 150.0 else 'Remission')
    
    flares = df_ibd[df_ibd['Flare_Status'] == 'Active_Flare']
    remissions = df_ibd[df_ibd['Flare_Status'] == 'Remission']
    
    print(f"   - IBD Samples with Calprotectin: Active Flare (>=150 ug/g) = {flares.shape[0]}, Remission (<150 ug/g) = {remissions.shape[0]}")
    
    flare_results = []
    
    if flares.shape[0] >= 3 and remissions.shape[0] >= 3:
        for col in species_cols:
            x = flares[col].values
            y = remissions[col].values
            
            u_stat, p_val = mannwhitneyu(x, y, alternative='two-sided')
            delta = calculate_cliffs_delta(x, y)
            
            mean_flare = np.mean(x)
            mean_remission = np.mean(y)
            
            flare_results.append({
                'Species': col.replace('_CLR', ''),
                'Mean_Active_Flare': mean_flare,
                'Mean_Remission': mean_remission,
                'Mean_Difference': mean_flare - mean_remission,
                'U_Statistic': u_stat,
                'p_value': p_val,
                'Cliffs_Delta': delta
            })
            
        flare_df = pd.DataFrame(flare_results)
        _, fdr_p, _, _ = multipletests(flare_df['p_value'], alpha=0.05, method='fdr_bh')
        flare_df['FDR_adjusted_p'] = fdr_p
        
        flare_df = flare_df.sort_values(by='FDR_adjusted_p')
        flare_output = os.path.join(RESULTS_DIR, 'flare_vs_remission_results.csv')
        flare_df.to_csv(flare_output, index=False)
        print(f"   - Saved Flare vs. Remission results to: {flare_output}")
    else:
        print("   - Warning: Insufficient sample size in either active flare or remission subgroups to run statistical tests.")
        flare_df = pd.DataFrame()

    # ==========================================
    # 🔍 4. Generating Publication-Quality Summary
    # ==========================================
    print("\n------------------------------------------------------------")
    print("Top Statistical Discoveries Summary:")
    print("------------------------------------------------------------")
    
    # Top Crohn's vs Healthy dysbiosis
    cd_sig = final_diff_df[(final_diff_df['Comparison'] == 'CD_vs_Healthy') & (final_diff_df['FDR_adjusted_p'] < 0.05)]
    print(f"\n1. Crohn's Disease Dysbiosis: {len(cd_sig)} species significantly altered (FDR < 0.05).")
    if not cd_sig.empty:
        # Show top 5 by effect size (Cliff's Delta)
        cd_sig_sorted = cd_sig.reindex(cd_sig['Cliffs_Delta'].abs().sort_values(ascending=False).index)
        for _, row in cd_sig_sorted.head(5).iterrows():
            dir_str = "depleted" if row['Cliffs_Delta'] < 0 else "enriched"
            print(f"   * {row['Species']}: {dir_str} (FDR p={row['FDR_adjusted_p']:.4e}, Cliff's Delta={row['Cliffs_Delta']:.3f})")

    # Top Ulcerative Colitis vs Healthy dysbiosis
    uc_sig = final_diff_df[(final_diff_df['Comparison'] == 'UC_vs_Healthy') & (final_diff_df['FDR_adjusted_p'] < 0.05)]
    print(f"\n2. Ulcerative Colitis Dysbiosis: {len(uc_sig)} species significantly altered (FDR < 0.05).")
    if not uc_sig.empty:
        uc_sig_sorted = uc_sig.reindex(uc_sig['Cliffs_Delta'].abs().sort_values(ascending=False).index)
        for _, row in uc_sig_sorted.head(5).iterrows():
            dir_str = "depleted" if row['Cliffs_Delta'] < 0 else "enriched"
            print(f"   * {row['Species']}: {dir_str} (FDR p={row['FDR_adjusted_p']:.4e}, Cliff's Delta={row['Cliffs_Delta']:.3f})")

    # Top Species correlated with Fecal Calprotectin
    corr_sig = corr_df[corr_df['FDR_adjusted_p'] < 0.05]
    print(f"\n3. Gut Inflammation (Calprotectin) Correlations: {len(corr_sig)} species significantly correlated (FDR < 0.05).")
    if not corr_sig.empty:
        for _, row in corr_sig.head(5).iterrows():
            rel_str = "positively correlated (marker of inflammation)" if row['Spearman_Rho'] > 0 else "negatively correlated (potentially protective/beneficial)"
            print(f"   * {row['Species']}: {rel_str} (Rho={row['Spearman_Rho']:.3f}, FDR p={row['FDR_adjusted_p']:.4e})")

    # Top Flare vs Remission biomarkers
    if not flare_df.empty:
        flare_sig = flare_df[flare_df['FDR_adjusted_p'] < 0.05]
        print(f"\n4. Active Flares vs. Remission Biomarkers: {len(flare_sig)} species differentially abundant (FDR < 0.05).")
        if not flare_sig.empty:
            for _, row in flare_sig.head(5).iterrows():
                dir_str = "depleted during flares" if row['Cliffs_Delta'] < 0 else "enriched during flares"
                print(f"   * {row['Species']}: {dir_str} (FDR p={row['FDR_adjusted_p']:.4e}, Cliff's Delta={row['Cliffs_Delta']:.3f})")
    print("------------------------------------------------------------")

if __name__ == "__main__":
    run_analysis()
