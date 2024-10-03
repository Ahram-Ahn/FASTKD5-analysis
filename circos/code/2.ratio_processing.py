import os
import pandas as pd
import numpy as np

samples = ['WT', 'KO']

gene_sequence = [
    "NCR1_TF",
    "RNR1",
    "TV",
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

# Genes to set values to zero, only mRNAs will be analyzed
genes_to_zero = [
    "RNR1",
    "RNR2",
    "ND6",
    "NCR",
    "NCR1",
    "NCR2",
    "un1",
    "un2",
    "un3",
    "un4"
    # 'tRNA' genes will be handled separately
]

# Path to the karyotype file (adjust if different)
karyotype_path = 'path/to/karyotype.txt'

# Base directories for WT and KO samples
base_graph_dirs = {
    'WT': "path/to/graph",
    'KO': "path/to/graph"
}

# File suffix for input CSV files
file_suffix = 'profile_masked_m-ratio-q0.csv'


def merge_gene_data(sample, base_graph_dir, gene_sequence, file_suffix):
    """
    Merges gene CSV files into consolidated plus and minus strand text files.
    Inserts the TV gene with zero counts after RNR1.
    """
    print(f"\nMerging gene data for sample: {sample}")
    
    # Define output file paths
    plus_output_path = os.path.join(base_graph_dir, f"{sample}_plus_ratio_merged.txt")
    minus_output_path = os.path.join(base_graph_dir, f"{sample}_minus_ratio_merged.txt")
    
    # tRNA-Val was missing in the data. To set the full length mt-genome, should be manually added with nan value 
    inserted_tv = False

    with open(plus_output_path, 'w') as plus_outfile, open(minus_output_path, 'w') as minus_outfile:
        for gene in gene_sequence:
            # Process Plus Strand 
            gene_plus_dir = os.path.join(base_graph_dir, gene)
            csv_plus_path = os.path.join(gene_plus_dir, 'full', file_suffix)
            
            if os.path.isfile(csv_plus_path):
                try:
                    df_plus = pd.read_csv(csv_plus_path)
                    # Check for required columns
                    if {'Position', 'Mutated'}.issubset(df_plus.columns):
                        for _, row in df_plus.iterrows():
                            position = row['Position']
                            count = row['Mutated']
                            plus_outfile.write(f"{gene} {position} {position} {count}\n")
                    else:
                        print(f"Warning: Required columns not found in {csv_plus_path}. Skipping plus strand for {gene}.")
                except Exception as e:
                    print(f"Error reading {csv_plus_path}: {e}. Skipping plus strand for {gene}.")
            else:
                print(f"Warning: {csv_plus_path} does not exist. Skipping plus strand for {gene}.")
            
            # Process Minus Strand
            gene_minus_dir = f"{gene}-minus"
            gene_minus_path = os.path.join(base_graph_dir, gene_minus_dir)
            csv_minus_path = os.path.join(gene_minus_path, 'full', file_suffix)
            
            if os.path.isfile(csv_minus_path):
                try:
                    df_minus = pd.read_csv(csv_minus_path)
                    # Check for required columns
                    if {'Position', 'Mutated'}.issubset(df_minus.columns):
                        min_position = df_minus['Position'].min()
                        for _, row in df_minus.iterrows():
                            position = row['Position']
                            adjusted_position = position - min_position + 1  # Adjust positions to start from 1
                            count = row['Mutated']
                            minus_outfile.write(f"{gene} {adjusted_position} {adjusted_position} {count}\n")
                    else:
                        print(f"Warning: Required columns not found in {csv_minus_path}. Skipping minus strand for {gene}.")
                except Exception as e:
                    print(f"Error reading {csv_minus_path}: {e}. Skipping minus strand for {gene}.")
            else:
                print(f"Warning: {csv_minus_path} does not exist. Skipping minus strand for {gene}.")
            
            # Insert TV Gene After RNR1
            if gene == "RNR1" and not inserted_tv:
                # Plus Strand
                for i in range(1, 71):  # 70 nucleotides
                    plus_outfile.write(f"TV {i} {i} 0\n")
                plus_outfile.write("\n\n")

                # Minus Strand
                for i in range(1, 71):
                    minus_outfile.write(f"TV {i} {i} 0\n")
                minus_outfile.write("\n\n")
                
                inserted_tv = True
                print("Inserted TV gene after RNR1.")
            
            if gene != "RNR1": 
                plus_outfile.write("\n\n")
                minus_outfile.write("\n\n")
    
    print(f"Merged data written to:\nPlus Strand: {plus_output_path}\nMinus Strand: {minus_output_path}")
    return plus_output_path, minus_output_path

def adjust_numbering(merged_path, adjusted_path):
    """
    Adjusts the Start and End numbering to start from 1.
    """
    print(f"Adjusting numbering for: {merged_path}")
    try:
        df = pd.read_csv(merged_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])
        df['Start'] = range(1, len(df) + 1)
        df['End'] = df['Start']
        df.to_csv(adjusted_path, sep="\t", index=False, header=False)
        print(f"Numbering adjusted and saved to: {adjusted_path}")
    except Exception as e:
        print(f"Error adjusting numbering for {merged_path}: {e}")

def modify_data(adjusted_plus_path, adjusted_minus_path, modified_plus_path, modified_minus_path, karyotype_df, genes_to_zero):
    """
    Maps positions to gene names, sets specified genes' values to zero, adds 0.001 to avoid division by zero, and rounds values.
    """
    print(f"\nModifying data for:\nPlus: {adjusted_plus_path}\nMinus: {adjusted_minus_path}")
    
    for strand, input_path, output_path in [
        ('Plus', adjusted_plus_path, modified_plus_path),
        ('Minus', adjusted_minus_path, modified_minus_path)
    ]:
        try:
            df = pd.read_csv(input_path, sep='\t', header=None, names=['Gene', 'Start', 'End', 'Value'])
            
            # Map positions to gene names
            df['Mapped_Gene'] = df['Start'].apply(lambda x: map_position_to_name(x, karyotype_df))
            
            # Replace original gene names with mapped ones
            df['Gene'] = df['Mapped_Gene']
            
            # Drop the temporary 'Mapped_Gene' column
            df.drop(columns=['Mapped_Gene'], inplace=True)
            
            # Set values to zero for specified genes
            df.loc[df['Gene'].isin(genes_to_zero), 'Value'] = 0.0
            
            # Identify and set tRNA genes to zero
            tRNA_genes = karyotype_df[karyotype_df['name'].str.startswith('tRNA')]['name'].unique()
            df.loc[df['Gene'].isin(tRNA_genes), 'Value'] = 0.0
            
            # Add 0.001 to all values to avoid division by zero in future calculations
            df['Value'] += 0.001
            
            # Round the 'Value' column to 3 decimal points
            df['Value'] = df['Value'].round(3)
            
            # Save the modified dataframe
            df.to_csv(output_path, sep="\t", index=False, header=False)
            print(f"Modified {strand} data saved to: {output_path}")
        
        except Exception as e:
            print(f"Error modifying {strand} data ({input_path}): {e}")

def map_position_to_name(pos, karyotype_df):
    """
    Maps a given position to a gene name using the karyotype dataframe.
    """
    matched = karyotype_df[(karyotype_df['start'] <= pos) & (pos <= karyotype_df['end'])]
    if not matched.empty:
        return matched.iloc[0]['name']
    else:
        return 'Unknown'

def compute_log_fold_change(wt_modified_plus, ko_modified_plus, final_output_path):
    """
    Computes the log2 fold change between KO and WT plus data and saves the result.
    """
    print(f"\nComputing log2 fold change between KO and WT.")
    
    try:
        # Read the modified plus data
        ko_plus_df = pd.read_csv(ko_modified_plus, sep='\t', header=None, names=['Gene', 'Start', 'End', 'KO_Value'])
        wt_plus_df = pd.read_csv(wt_modified_plus, sep='\t', header=None, names=['Gene', 'Start', 'End', 'WT_Value'])
        
        # Merge the data on 'Gene', 'Start', and 'End'
        merged_plus_df = pd.merge(ko_plus_df, wt_plus_df, on=['Gene', 'Start', 'End'])
        
        # Compute log2 fold difference
        merged_plus_df['Log2_Fold_Change'] = np.log2(merged_plus_df['KO_Value'] / merged_plus_df['WT_Value'])
        
        # Handle possible infinite or NaN values
        merged_plus_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        merged_plus_df.dropna(subset=['Log2_Fold_Change'], inplace=True)
        
        # Round the log fold change to three decimal points
        merged_plus_df['Log2_Fold_Change'] = merged_plus_df['Log2_Fold_Change'].round(3)
        
        # Save the final data
        merged_plus_df[['Gene', 'Start', 'End', 'Log2_Fold_Change']].to_csv(final_output_path, sep='\t', index=False, header=False)
        
        print(f"Log2 fold difference computed and saved to: {final_output_path}")
    
    except Exception as e:
        print(f"Error computing log2 fold change: {e}")

def load_karyotype(karyotype_path):
    """
    Loads the karyotype data into a dataframe.
    """
    try:
        karyotype_cols = ['chr', 'dash', 'name', 'abbrev', 'start', 'end', 'chromosome']
        karyotype = pd.read_csv(karyotype_path, sep='\s+', header=None, names=karyotype_cols)
        print(f"Karyotype data loaded from: {karyotype_path}")
        return karyotype
    except Exception as e:
        print(f"Error loading karyotype data: {e}")
        return None

def main():
    # Load karyotype data once
    karyotype_df = load_karyotype(karyotype_path)
    if karyotype_df is None:
        print("Karyotype data is essential for processing. Exiting.")
        return
    
    # Dictionary to store paths for each sample
    processed_paths = {}
    
    for sample in samples:
        print(f"\n===== Processing Sample: {sample} =====")
        base_graph_dir = base_graph_dirs[sample]
        
        # Merge Gene Data
        merged_plus, merged_minus = merge_gene_data(sample, base_graph_dir, gene_sequence, file_suffix)
        
        # Adjust Numbering
        adjusted_plus = os.path.join(base_graph_dir, f"{sample}_plus_ratio_adjusted_numbering.txt")
        adjusted_minus = os.path.join(base_graph_dir, f"{sample}_minus_ratio_adjusted_numbering.txt")
        adjust_numbering(merged_plus, adjusted_plus)
        adjust_numbering(merged_minus, adjusted_minus)
        
        # Modify Data
        modified_plus = os.path.join(base_graph_dir, f"{sample}_plus_data_modified_final.txt")
        modified_minus = os.path.join(base_graph_dir, f"{sample}_minus_data_modified_final.txt")
        modify_data(adjusted_plus, adjusted_minus, modified_plus, modified_minus, karyotype_df, genes_to_zero)
        
        # Store modified paths for later use
        processed_paths[sample] = {
            'plus': modified_plus,
            'minus': modified_minus
        }
    
    # Compute Log2 Fold Change (KO vs WT)
    if 'WT' in processed_paths and 'KO' in processed_paths:
        final_plus_output = os.path.join(
            base_graph_dirs['KO'], "WT_vs_KO_log_fold_difference_plus_final.txt"
        )
        compute_log_fold_change(
            wt_modified_plus=processed_paths['WT']['plus'],
            ko_modified_plus=processed_paths['KO']['plus'],
            final_output_path=final_plus_output
        )
    else:
        print("Both WT and KO samples must be processed to compute log2 fold changes.")

if __name__ == "__main__":
    main()
