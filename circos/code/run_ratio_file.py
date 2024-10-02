import os
import pandas as pd

# Define the updated karyotype as a multi-line string
karyotype = """
chr - NCR1 NCR 0 576 lgrey
chr - tRNAF F 577 647 chr4
chr - RNR1 RNR1 648 1601 chr8
chr - tRNAV V 1602 1670 chr4
chr - RNR2 RNR2 1671 3229 chr8
chr - tRNAL1 L1 3230 3305 chr4
chr - ND1 ND1 3306 4262 chr15
chr - tRNAI I 4263 4331 chr4
chr - tRNAQ Q 4332 4401 chr4
chr - tRNAM M 4402 4469 chr4
chr - ND2 ND2 4470 5511 chr15
chr - tRNAW W 5512 5579 chr4
chr - tRNAA A 5587 5656 chr4
chr - tRNAN N 5657 5729 chr4
chr - un1 - 5730 5760 white
chr - tRNAC C 5761 5826 chr4
chr - tRNAY Y 5827 5891 chr4
chr - COX1 COX1 5902 7445 chr15
chr - tRNAS1 S1 7446 7517 chr4
chr - tRNAD D 7518 7585 chr4
chr - COX2 COX2 7586 8269 chr15
chr - un2 - 8270 8294 white
chr - tRNAK K 8295 8365 chr4
chr - tRNAATP8/6 ATP8/6 8366 9207 chr15
chr - COX3 COX3 9208 9990 chr15
chr - tRNAG G 9991 10058 chr4
chr - ND3 ND3 10059 10404 chr15
chr - tRNAR R 10405 10469 chr4
chr - ND4L/4 ND4L/4 10470 12137 chr15
chr - tRNAH H 12138 12206 chr4
chr - tRNAS2 S2 12207 12265 chr4
chr - tRNAL2 L2 12266 12336 chr4
chr - ND5 ND5 12337 14148 chr15
chr - ND6 ND6 14149 14673 chr15
chr - tRNAE E 14674 14746 chr4
chr - CYTB CYTB 14747 15887 chr15
chr - tRNAT T 15888 15955 chr4
chr - tRNAP P 15956 16023 chr4
chr - NCR2 NCR 16024 16569 lgrey
"""
import os
import pandas as pd

def main():
    # Define the base graph directory. beaware of the directory
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
    plus_output_path = os.path.join(base_graph_dir, "KO_plus_ratio_merged.txt")
    minus_output_path = os.path.join(base_graph_dir, "KO_minus_ratio_merged.txt")
    
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
                        for _, row in df_minus.iterrows():
                            position = row['Position']
                            count = row['Mutated']
                            # Negate the count for minus strand
                            adjusted_count = count
                            # Write: gene, position, position, adjusted_count
                            minus_outfile.write(f"{gene} {position} {position} {adjusted_count}\n")
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
