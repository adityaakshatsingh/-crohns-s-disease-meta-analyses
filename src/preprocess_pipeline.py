"""
Preprocessing Pipeline for HMP2 IBD Metagenomics Dataset
=========================================================

Author: Aditya Akshat Singh (CSE Student) & Antigravity (AI Coding Assistant)
Date: June 2026
Description: This script performs clean deduplication, chronological sorting, 
             rare species prevalence filtering, participant criteria filtering,
             biweekly time-series grid resampling/interpolation, and 
             compositional Centered Log-Ratio (CLR) transformation.
             
No external dependencies (runs in vanilla Python with default libraries).
"""

import csv
import os
import math
from collections import defaultdict

# Define file paths
CLEANED_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_cleaned.csv"
OUTPUT_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_preprocessed.csv"

def run_preprocessing():
    print("--------------------------------------------------")
    print("Starting Metagenomics Preprocessing Pipeline...")
    print("--------------------------------------------------")

    if not os.path.exists(CLEANED_FILE):
        print(f"Cleaned dataset not found at: {CLEANED_FILE}")
        print("Please run clean_dataset.py first.")
        return

    # Read cleaned data
    with open(CLEANED_FILE, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    num_samples = len(rows)
    species_cols = header[5:]
    num_species = len(species_cols)

    # ==========================================
    # Step 1: Prevalence-based Species Filtering
    # ==========================================
    print("Step 1: Filtering rare microbial species...")
    kept_species_indices = []
    for idx in range(5, len(header)):
        col_name = header[idx]
        values = [float(r[idx]) if r[idx] else 0.0 for r in rows]
        
        prevalence = sum(1 for v in values if v > 0.0) / num_samples
        mean_abundance = sum(values) / num_samples
        
        # Prevalence > 10% AND mean abundance > 0.01%
        if prevalence > 0.10 and mean_abundance > 0.01:
            kept_species_indices.append(idx)

    print(f"   - Kept {len(kept_species_indices)} out of {num_species} species.")

    # ==========================================
    # Step 2: Participant Criteria Filtering
    # ==========================================
    print("Step 2: Filtering participants by clinical history coverage...")
    participants_data = defaultdict(list)
    for r in rows:
        pid = r[1]
        participants_data[pid].append(r)

    kept_participants = {}
    for pid, p_rows in participants_data.items():
        total = len(p_rows)
        missing_fcal = sum(1 for r in p_rows if not r[4])
        missing_pct = missing_fcal / total if total > 0 else 1.0
        
        # Criteria: Calprotectin missingness <= 70% and total samples >= 4
        if missing_pct <= 0.70 and total >= 4:
            kept_participants[pid] = p_rows

    print(f"   - Kept {len(kept_participants)} out of {len(participants_data)} participants.")

    # Helper function for linear interpolation / nearest neighbor extrapolation
    def interpolate_val(w, time_val_pairs):
        valid_pairs = [(t, v) for t, v in time_val_pairs if v is not None]
        if not valid_pairs:
            return None
        
        # Exact match
        for t, v in valid_pairs:
            if abs(t - w) < 1e-6:
                return v
                
        # Find closest indices before and after
        before = [(t, v) for t, v in valid_pairs if t < w]
        after = [(t, v) for t, v in valid_pairs if t > w]
        
        if before and after:
            t1, v1 = before[-1]
            t2, v2 = after[0]
            return v1 + (w - t1) / (t2 - t1) * (v2 - v1)
        elif before:
            return before[-1][1]  # Extrapolate nearest neighbor (carry forward)
        elif after:
            return after[0][1]   # Extrapolate nearest neighbor (carry backward)
        return None

    # ==========================================================
    # Step 3 & 4: Resampling, Interpolation & CLR Transformation
    # ==========================================================
    print("Step 3: Resampling onto biweekly grid & interpolating...")
    print("Step 4: Applying Centered Log-Ratio (CLR) transformation...")
    
    preprocessed_rows = []

    for pid, p_rows in kept_participants.items():
        # Sort chronologically
        p_rows.sort(key=lambda x: float(x[2]) if x[2] else 0.0)
        
        diagnosis = p_rows[0][3]
        times = [float(r[2]) if r[2] else 0.0 for r in p_rows]
        min_week = min(times)
        max_week = max(times)
        
        # Target grid (even integer weeks)
        start_grid = int(math.ceil(min_week / 2.0) * 2)
        end_grid = int(math.floor(max_week / 2.0) * 2)
        grid = list(range(start_grid, end_grid + 1, 2))
        
        # Compile time-value pairs for Calprotectin
        fcal_pairs = []
        for r in p_rows:
            t = float(r[2]) if r[2] else 0.0
            v = float(r[4]) if r[4] else None
            fcal_pairs.append((t, v))
        fcal_pairs.sort(key=lambda x: x[0])
        
        # Compile time-value pairs for species
        species_pairs = defaultdict(list)
        for r in p_rows:
            t = float(r[2]) if r[2] else 0.0
            for s_idx in kept_species_indices:
                v = float(r[s_idx]) if r[s_idx] else 0.0
                species_pairs[s_idx].append((t, v))
                
        for s_idx in kept_species_indices:
            species_pairs[s_idx].sort(key=lambda x: x[0])
            
        # Interpolate and transform each grid week
        for w in grid:
            w_float = float(w)
            
            # 1. Calprotectin
            interp_fcal = interpolate_val(w_float, fcal_pairs)
            fcal_str = f"{interp_fcal:.4f}" if interp_fcal is not None else ""
            
            # 2. Species abundances
            abundances = []
            for s_idx in kept_species_indices:
                v_interp = interpolate_val(w_float, species_pairs[s_idx])
                abundances.append(v_interp if v_interp is not None else 0.0)
                
            # 3. CLR transformation
            # Add a 1e-6 pseudocount to break zero composition limits
            pseudocount = 1e-6
            abundances_pseudo = [a + pseudocount for a in abundances]
            
            logs = [math.log(a) for a in abundances_pseudo]
            mean_log = sum(logs) / len(logs)
            clr_values = [l - mean_log for l in logs]
            
            # Format output decimals
            clr_str_values = [f"{val:.6f}" for val in clr_values]
            
            # Assemble preprocessed row
            ext_id = f"{pid}_W{w}"
            row_out = [ext_id, pid, f"{w_float:.1f}", diagnosis, fcal_str] + clr_str_values
            preprocessed_rows.append(row_out)

    # ==========================
    # Step 5: Save Preprocessed Data
    # ==========================
    print("Step 5: Saving preprocessed dataset...")
    out_header = ["External ID", "Participant ID", "week_num", "diagnosis", "fecalcal"]
    for s_idx in kept_species_indices:
        out_header.append(f"{header[s_idx]}_CLR")

    # Final sort for easy viewing
    preprocessed_rows.sort(key=lambda x: (x[1], float(x[2])))

    with open(OUTPUT_FILE, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(out_header)
        writer.writerows(preprocessed_rows)

    print("--------------------------------------------------")
    print("Preprocessing Pipeline Completed Successfully!")
    print(f"   - Final sample rows: {len(preprocessed_rows)}")
    print(f"   - Final columns: {len(out_header)}")
    print(f"   - Output saved to: {OUTPUT_FILE}")
    print("--------------------------------------------------")

if __name__ == "__main__":
    run_preprocessing()
