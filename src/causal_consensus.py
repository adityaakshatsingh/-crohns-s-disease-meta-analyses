"""
Consensus Causal Network Construction & Validation
===================================================

Author: Aditya Akshat Singh (CSE Student) & Antigravity (AI Coding Assistant)
Date: June 2026
Description: This script executes the final Phase 6 of the research plan:
             1. Load causal results from Granger, DBN, and PCMCI (Tigramite).
             2. Matches directed edges across the three independent methods.
             3. Generates a consensus dataset showing the intersection of causal support.
             4. Outputs the final consensus table to results/causal_consensus_edges.csv.

Dependencies: pandas, numpy
"""

import os
import pandas as pd

# File paths
GRANGER_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\results\causal_granger_results.csv"
DBN_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\results\causal_dbn_results.csv"
TIGRAMITE_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\results\causal_tigramite_results.csv"
OUTPUT_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\results\causal_consensus_edges.csv"

def run_consensus_analysis():
    print("------------------------------------------------------------")
    print("Starting Phase 6: Consensus Causal Network Analysis...")
    print("------------------------------------------------------------")

    # Verify input files exist
    files = [GRANGER_FILE, DBN_FILE, TIGRAMITE_FILE]
    for f in files:
        if not os.path.exists(f):
            print(f"Error: Required file not found: {f}")
            return

    # Load results
    granger_df = pd.read_csv(GRANGER_FILE)
    dbn_df = pd.read_csv(DBN_FILE)
    tig_df = pd.read_csv(TIGRAMITE_FILE)

    consensus_dict = {}

    # Helper function to add/update edges in consensus
    def add_edge(pred, target, method, p_val, weight):
        key = (pred, target)
        if key not in consensus_dict:
            consensus_dict[key] = {
                'Predictor': pred,
                'Target': target,
                'Granger_Support': 0,
                'Granger_p': None,
                'DBN_Support': 0,
                'DBN_p': None,
                'DBN_Coef': None,
                'PCMCI_Support': 0,
                'PCMCI_p': None,
                'PCMCI_MCI': None,
                'Consensus_Score': 0
            }
        
        c = consensus_dict[key]
        if method == 'Granger':
            c['Granger_Support'] = 1
            c['Granger_p'] = p_val
        elif method == 'DBN':
            c['DBN_Support'] = 1
            c['DBN_p'] = p_val
            c['DBN_Coef'] = weight
        elif method == 'PCMCI':
            c['PCMCI_Support'] = 1
            c['PCMCI_p'] = p_val
            c['PCMCI_MCI'] = weight

    # 1. Parse Granger results (FDR < 0.05)
    for _, row in granger_df.iterrows():
        sp = row['Species']
        # Forward: Species -> Calprotectin
        if row['Forward_FDR_p'] < 0.05:
            add_edge(sp, 'fecalcal', 'Granger', row['Forward_FDR_p'], None)
        # Backward: Calprotectin -> Species
        if row['Backward_FDR_p'] < 0.05:
            add_edge('fecalcal', sp, 'Granger', row['Backward_FDR_p'], None)

    # 2. Parse DBN results (FDR < 0.05, exclude self-loops)
    for _, row in dbn_df.iterrows():
        pred = row['Predictor']
        tar = row['Target']
        if pred == tar:
            continue
        if row['FDR_adjusted_p'] < 0.05:
            add_edge(pred, tar, 'DBN', row['FDR_adjusted_p'], row['Coefficient'])

    # 3. Parse PCMCI/Tigramite results (FDR < 0.05, exclude self-loops)
    for _, row in tig_df.iterrows():
        pred = row['Predictor']
        tar = row['Target']
        if pred == tar:
            continue
        if row['FDR_adjusted_p'] < 0.05:
            add_edge(pred, tar, 'PCMCI', row['FDR_adjusted_p'], row['MCI_Correlation'])

    # 4. Calculate Consensus Scores
    consensus_list = list(consensus_dict.values())
    for c in consensus_list:
        score = c['Granger_Support'] + c['DBN_Support'] + c['PCMCI_Support']
        c['Consensus_Score'] = score

    # Convert to DataFrame
    consensus_df = pd.DataFrame(consensus_list)
    consensus_df = consensus_df.sort_values(by='Consensus_Score', ascending=False)

    # Save to CSV
    consensus_df.to_csv(OUTPUT_FILE, index=False)
    print(f"   - Successfully saved consensus causal edges to: {OUTPUT_FILE}")

    # Print summary
    print("\n------------------------------------------------------------")
    print("Consensus Directed Causal Edges (Score >= 2):")
    print("------------------------------------------------------------")
    score_3 = consensus_df[consensus_df['Consensus_Score'] == 3]
    score_2 = consensus_df[consensus_df['Consensus_Score'] == 2]

    print(f"Total 3-Method Consensus Edges: {len(score_3)}")
    for _, row in score_3.iterrows():
        rel_str = "promotes" if (row['DBN_Coef'] and row['DBN_Coef'] > 0) or (row['PCMCI_MCI'] and row['PCMCI_MCI'] > 0) else "inhibits"
        print(f"   * {row['Predictor']} ---> {row['Target']} [3-Method Consensus] ({rel_str})")
        print(f"     - Granger p={row['Granger_p']:.2e} | DBN Coef={row['DBN_Coef']:.3f} | PCMCI MCI={row['PCMCI_MCI']:.3f}")

    print(f"\nTotal 2-Method Consensus Edges: {len(score_2)}")
    for _, row in score_2.iterrows():
        # Determine relationship direction string
        coef_val = row['DBN_Coef'] if pd.notna(row['DBN_Coef']) else row['PCMCI_MCI']
        rel_str = "promotes" if (coef_val and coef_val > 0) else "inhibits"
        
        methods = []
        if row['Granger_Support'] == 1: methods.append('Granger')
        if row['DBN_Support'] == 1: methods.append('DBN')
        if row['PCMCI_Support'] == 1: methods.append('PCMCI')
        
        print(f"   * {row['Predictor']} ---> {row['Target']} [2-Method Consensus: {', '.join(methods)}] ({rel_str})")
    print("------------------------------------------------------------")

if __name__ == "__main__":
    run_consensus_analysis()
