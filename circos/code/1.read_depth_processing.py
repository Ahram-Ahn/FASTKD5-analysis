import os
import pandas as pd

def merge_gene_data(base_graph_dir, gene_sequence, output_prefix, data_type='Unambiguous', file_suffix='profile_related_n-count.csv'):
    """
    Merges data from multiple gene files into plus and minus strand files.

    Parameters:
    - base_graph_dir: Base directory containing gene data.
    - gene_sequence: List of gene names.
    - output_prefix: Prefix for output files.
    - data_type: Column name to extract from gene data files.
    - file_suffix: Suffix of the gene data files.
    """
    plus_output_path = os.path.join(base_graph_dir, f"{output_prefix}_plus_merged.txt")
    minus_output_path = os.path.join(base_graph_dir, f"{output_prefix}_minus_merged.txt")

    # TV (tRNA-Val) should be manually added to the processed data
    # NCR1_TF did not pass the siesmic-rna in the minus strand, so should be manually added
 
    inserted_tv_plus = False
    inserted_tv_minus = False
    ncr1_tf_length_plus = 0

    with open(plus_output_path, 'w') as plus_outfile, open(minus_output_path, 'w') as minus_outfile:
        for gene in gene_sequence:
            # Process Plus Strand
            gene_plus_dir = os.path.join(base_graph_dir, gene)
            csv_plus_path = os.path.join(gene_plus_dir, 'full', file_suffix)

            if os.path.isfile(csv_plus_path):
                try:
                    df_plus = pd.read_csv(csv_plus_path)
                    if {'Position', data_type}.issubset(df_plus.columns):
                        for _, row in df_plus.iterrows():
                            position = row['Position']
                            count = row[data_type]
                            plus_outfile.write(f"{gene} {position} {position} {count}\n")
                        # If gene is NCR1_TF, store its length for minus strand
                        if gene == "NCR1_TF":
                            ncr1_tf_length_plus = df_plus['Position'].max() 
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

            if gene == "NCR1_TF":
                # Special handling for NCR1_TF in minus strand
                if os.path.isfile(csv_minus_path):
                    try:
                        df_minus = pd.read_csv(csv_minus_path)
                        if {'Position', data_type}.issubset(df_minus.columns):
                            for _, row in df_minus.iterrows():
                                position = row['Position']
                                adjusted_position = position - df_minus['Position'].min() + 1  # To make sure positions start from 1
                                count = row[data_type]
                                minus_outfile.write(f"{gene} {adjusted_position} {adjusted_position} {count}\n")
                        else:
                            print(f"Warning: Required columns not found in {csv_minus_path}. Inserting NCR1_TF with zero counts in minus strand.")
                            # Insert NCR1_TF with zero counts based on plus strand length
                            if ncr1_tf_length_plus > 0:
                                for i in range(1, ncr1_tf_length_plus + 1):
                                    minus_outfile.write(f"{gene} {i} {i} 0\n")
                            else:
                                print("Error: Cannot determine NCR1_TF length from plus strand. Skipping insertion.")
                    except Exception as e:
                        print(f"Error reading {csv_minus_path}: {e}. Inserting NCR1_TF with zero counts in minus strand.")
                        # Insert NCR1_TF with zero counts based on plus strand length
                        if ncr1_tf_length_plus > 0:
                            for i in range(1, ncr1_tf_length_plus + 1):
                                minus_outfile.write(f"{gene} {i} {i} 0\n")
                        else:
                            print("Error: Cannot determine NCR1_TF length from plus strand. Skipping insertion.")
                else:
                    print(f"Warning: {csv_minus_path} does not exist. Inserting NCR1_TF with zero counts in minus strand.")
                    # Insert NCR1_TF with zero counts based on plus strand length
                    if ncr1_tf_length_plus > 0:
                        for i in range(1, ncr1_tf_length_plus + 1):
                            minus_outfile.write(f"{gene} {i} {i} 0\n")
                    else:
                        print("Error: Cannot determine NCR1_TF length from plus strand. Skipping insertion.")
            else:
                # Regular processing for other genes
                if os.path.isfile(csv_minus_path):
                    try:
                        df_minus = pd.read_csv(csv_minus_path)
                        if {'Position', data_type}.issubset(df_minus.columns):
                            min_position = df_minus['Position'].min()
                            for _, row in df_minus.iterrows():
                                position = row['Position']
                                adjusted_position = position - min_position + 1  # Adjust positions to start from 1
                                count = row[data_type]
                                minus_outfile.write(f"{gene} {adjusted_position} {adjusted_position} {count}\n")
                        else:
                            print(f"Warning: Required columns not found in {csv_minus_path}. Skipping minus strand for {gene}.")
                    except Exception as e:
                        print(f"Error reading {csv_minus_path}: {e}. Skipping minus strand for {gene}.")
                else:
                    print(f"Warning: {csv_minus_path} does not exist. Skipping minus strand for {gene}.")

            # Insert TV after RNR1
            if gene == "RNR1" and not inserted_tv_plus and not inserted_tv_minus:
                # Insert TV in Plus Strand
                for i in range(1, 71):  # 70 nucleotides
                    plus_outfile.write(f"TV {i} {i} 0\n")
                plus_outfile.write("\n\n")  # Blank lines after gene's data
                inserted_tv_plus = True

                # Insert TV in Minus Strand
                for i in range(1, 71):
                    minus_outfile.write(f"TV {i} {i} 0\n")
                minus_outfile.write("\n\n")
                inserted_tv_minus = True

            # Insert Blank Lines After Each Gene's Data (excluding the special TV insertion)
            if gene != "RNR1":  # TV already added with blank lines
                plus_outfile.write("\n\n")
                minus_outfile.write("\n\n")

    print(f"Merged data has been written to:\nPlus Strand: {plus_output_path}\nMinus Strand: {minus_output_path}")
    return plus_output_path, minus_output_path

def normalize_data(plus_path, minus_path, output_prefix):
    """
    Normalizes the merged data based on total counts and adjusts numbering.

    Parameters:
    - plus_path: Path to the plus strand merged data file.
    - minus_path: Path to the minus strand merged data file.
    - output_prefix: Prefix for output files.
    """
    # Read data
    try:
        plus_df = pd.read_csv(plus_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])
    except Exception as e:
        print(f"Error reading plus strand file {plus_path}: {e}")
        return None, None

    try:
        minus_df = pd.read_csv(minus_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])
    except Exception as e:
        print(f"Error reading minus strand file {minus_path}: {e}")
        return None, None

    # Compute total sum
    total_sum = plus_df['Value'].sum() + minus_df['Value'].sum()

    if total_sum == 0:
        print("Total sum of counts is zero. Cannot normalize.")
        return None, None

    # Normalize
    plus_df['Value'] = round((plus_df['Value'] * 16569) / total_sum, 3)
    minus_df['Value'] = round((minus_df['Value'] * 16569) / total_sum, 3)

    # Adjust numbering to be from 1 to the length of the whole sequence
    plus_df['Start'] = range(1, len(plus_df) + 1)
    plus_df['End'] = plus_df['Start']
    minus_df['Start'] = range(1, len(minus_df) + 1)
    minus_df['End'] = minus_df['Start']

    # Save normalized data
    output_dir = 'path/to/output_dir'
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    plus_output = os.path.join(output_dir, f"{output_prefix}_plus_normalized.txt")
    minus_output = os.path.join(output_dir, f"{output_prefix}_minus_normalized.txt")

    try:
        plus_df.to_csv(plus_output, sep="\t", index=False, header=False)
        print(f"Normalized plus strand data saved to: {plus_output}")
    except Exception as e:
        print(f"Error saving normalized plus strand data: {e}")

    try:
        minus_df.to_csv(minus_output, sep="\t", index=False, header=False)
        print(f"Normalized minus strand data saved to: {minus_output}")
    except Exception as e:
        print(f"Error saving normalized minus strand data: {e}")

    return plus_output, minus_output

if __name__ == "__main__":
    # Define the base graph directory (replace with your own path)
    base_graph_dir = "/path/to/graph"

    # Define the gene sequence in the specified order (order given by the karyotype)
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

    # Process both WT and KO samples
    for sample in ['WT', 'KO']:
        print(f"\nProcessing sample: {sample}")
        # Merge data
        plus_merged_path, minus_merged_path = merge_gene_data(
            base_graph_dir, gene_sequence, output_prefix=sample, data_type='Unambiguous', file_suffix='profile_related_n-count.csv'
        )

        # Normalize data
        normalize_data(plus_merged_path, minus_merged_path, output_prefix=sample)
