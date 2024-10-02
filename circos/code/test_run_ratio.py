import os
import pandas as pd

def main():
    # Define the base graph directory. Be aware of the directory
    base_graph_dir = "/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph"
    
    # Define the gene sequence in the specified order
    gene_sequence = [
        "NCR1_TF",
        "RNR1",
        "RNR2",
        "TL1_ND1",
        "TI_TQ_ND2",
        "TW_TA_TN_TC_TY_COX1",
        "TD_COX2",
        "TK_ATP86_COX3",
        "TG_ND3",
        "TR_ND4L4_TH_TS2_TL2",
        "ND5_CYTB",
        "TT_TP_NCR2"
    ]
    
    # Initialize the output file paths
    plus_output_path = os.path.join(base_graph_dir, "KO_plus_ratio_merged_test.txt")
    minus_output_path = os.path.join(base_graph_dir, "KO_minus_ratio_merged_test.txt")
    
    # Open the output files for writing
    with open(plus_output_path, 'w') as plus_outfile, open(minus_output_path, 'w') as minus_outfile:
        # Iterate through each gene in the specified sequence
        for gene in gene_sequence:
            # --- Process Plus Strand ---
            gene_plus_dir = os.path.join(base_graph_dir, gene)
            csv_plus_path = os.path.join(gene_plus_dir, 'full', 'profile_masked_m-ratio-q0.csv')
            
            if os.path.isfile(csv_plus_path):
                try:
                    df_plus = pd.read_csv(csv_plus_path)
                    # Check for required columns
                    if {'Position', 'Mutated'}.issubset(df_plus.columns):
                        for _, row in df_plus.iterrows():
                            position = row['Position']
                            count = row['Mutated']
                            # Write: gene, position, position, count
                            plus_outfile.write(f"{gene} {position} {position} {count}\n")
                    else:
                        print(f"Warning: Required columns not found in {csv_plus_path}. Skipping plus strand for {gene}.")
                except Exception as e:
                    print(f"Error reading {csv_plus_path}: {e}. Skipping plus strand for {gene}.")
            else:
                print(f"Warning: {csv_plus_path} does not exist. Skipping plus strand for {gene}.")
            
            # --- Process Minus Strand ---
            gene_minus_dir = f"{gene}-minus"
            gene_minus_path = os.path.join(base_graph_dir, gene_minus_dir)
            csv_minus_path = os.path.join(gene_minus_path, 'full', 'profile_masked_m-ratio-q0.csv')
            
            if os.path.isfile(csv_minus_path):
                try:
                    df_minus = pd.read_csv(csv_minus_path)
                    # Check for required columns
                    if {'Position', 'Mutated'}.issubset(df_minus.columns):
                        # Adjust positions to start from 1
                        min_position = df_minus['Position'].min()
                        for _, row in df_minus.iterrows():
                            position = row['Position']
                            adjusted_position = position - min_position + 1  # Start from 1
                            count = row['Mutated']
                            # Write: gene, adjusted_position, adjusted_position, count
                            minus_outfile.write(f"{gene} {adjusted_position} {adjusted_position} {count}\n")
                    else:
                        print(f"Warning: Required columns not found in {csv_minus_path}. Skipping minus strand for {gene}.")
                except Exception as e:
                    print(f"Error reading {csv_minus_path}: {e}. Skipping minus strand for {gene}.")
            else:
                print(f"Warning: {csv_minus_path} does not exist. Skipping minus strand for {gene}.")
            
            # --- Insert Blank Lines After Each Gene's Data ---
            plus_outfile.write("\n\n")
            minus_outfile.write("\n\n")
    
    print(f"Merged data has been written to:")
    print(f"Plus Strand: {plus_output_path}")
    print(f"Minus Strand: {minus_output_path}")

if __name__ == "__main__":
    main()
