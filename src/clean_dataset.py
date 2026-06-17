"""
Dataset Cleaning and Deduplication for HMP2 IBD Metagenomics Dataset
===================================================================

Author: Aditya Akshat Singh (CSE Student) & Antigravity (AI Coding Assistant)
Date: June 2026
Description: This script cleans the raw HMP2 Metagenomics Atlas CSV file:
             1. Deduplicates exact matching rows (which make up ~60% of the dataset).
             2. Resolves sample duplication: groups by External ID (sample ID) 
                and retains the row that contains Fecal Calprotectin data.
             3. Sorts all records chronologically by Participant ID and week.
             4. Merges duplicate weeks for the same participant by averaging values.
"""

import csv
import os
from collections import defaultdict

# File paths
RAW_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_atlas_20260219_121629.csv"
OUTPUT_FILE = r"C:\Users\voltt\OneDrive\Desktop\heal\data\hmp2_ibd_metagenomics_cleaned.csv"

def run_cleaning():
    print("--------------------------------------------------")
    print("Starting Dataset Cleaning and Deduplication...")
    print("--------------------------------------------------")

    if not os.path.exists(RAW_FILE):
        print(f"Error: Raw file not found at: {RAW_FILE}")
        return

    # Read raw data
    with open(RAW_FILE, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    print(f"Read {len(rows)} total rows from raw file.")

    # Step 1: Group rows by External ID (Sample ID)
    samples_by_id = defaultdict(list)
    for r in rows:
        samples_by_id[r[0]].append(r)

    print(f"Found {len(samples_by_id)} unique sample IDs (External IDs).")

    unique_samples = []
    # Step 2: Keep the best row for each sample ID (prioritize non-empty calprotectin)
    for ext_id, group in samples_by_id.items():
        # Sort so rows with calprotectin (index 4 populated) come first
        group.sort(key=lambda x: 0 if x[4] else 1)
        unique_samples.append(group[0])

    # Step 3: Group by Participant ID for sorting and week resolution
    participants_data = defaultdict(list)
    for r in unique_samples:
        pid = r[1]
        participants_data[pid].append(r)

    cleaned_dataset = []

    # Step 4: Resolve multiple samples in the same week for a participant
    for pid, p_rows in participants_data.items():
        # Sort chronologically by week_num
        parsed_rows = []
        for r in p_rows:
            try:
                week = float(r[2]) if r[2] else 0.0
            except ValueError:
                week = 0.0
            parsed_rows.append((week, r))
        
        parsed_rows.sort(key=lambda x: x[0])
        
        # Group by week_num
        weeks_dict = defaultdict(list)
        for week, r in parsed_rows:
            weeks_dict[week].append(r)
            
        for week, week_rows in sorted(weeks_dict.items()):
            if len(week_rows) == 1:
                cleaned_dataset.append(week_rows[0])
            else:
                # Merge duplicate weeks by averaging abundances and calprotectin
                print(f"   - Merging duplicate week {week} for participant {pid} (Count: {len(week_rows)})")
                
                template = week_rows[0]
                ext_id = template[0]
                diagnosis = template[3]
                
                # Calprotectin values
                calprotectins = []
                for r in week_rows:
                    if r[4]:
                        try:
                            calprotectins.append(float(r[4]))
                        except ValueError:
                            pass
                avg_calprotectin = str(sum(calprotectins) / len(calprotectins)) if calprotectins else ""
                
                # Species abundances
                num_species = len(header) - 5
                abundances_sum = [0.0] * num_species
                for r in week_rows:
                    for idx in range(5, len(header)):
                        val = float(r[idx]) if r[idx] else 0.0
                        abundances_sum[idx - 5] += val
                
                avg_abundances = [str(x / len(week_rows)) for x in abundances_sum]
                
                merged_row = [ext_id, pid, str(week), diagnosis, avg_calprotectin] + avg_abundances
                cleaned_dataset.append(merged_row)

    # Step 5: Final sort by Participant ID and week_num
    def get_sort_key(row):
        try:
            week = float(row[2])
        except ValueError:
            week = 0.0
        return (row[1], week)

    cleaned_dataset.sort(key=get_sort_key)

    # Write cleaned data
    with open(OUTPUT_FILE, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(cleaned_dataset)

    print("--------------------------------------------------")
    print("Dataset Cleaning Completed Successfully!")
    print(f"   - Cleaned rows: {len(cleaned_dataset)}")
    print(f"   - Output saved to: {OUTPUT_FILE}")
    print("--------------------------------------------------")

if __name__ == "__main__":
    run_cleaning()
