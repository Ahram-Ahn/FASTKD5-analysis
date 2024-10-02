import pandas as pd

# Load the data
minus_data_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_minus_ratio_merged.txt'
plus_data_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_plus_ratio_merged.txt'

# Read the files into pandas dataframes
minus_df = pd.read_csv(minus_data_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])
plus_df = pd.read_csv(plus_data_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])

# Adjust numbering
def adjust_numbering(df, start_index):
    # Create contiguous numbering for 'Start' and 'End'
    df['Start'] = range(start_index, start_index + len(df))
    df['End'] = df['Start']
    return df, start_index + len(df)

# Adjust numbering for plus and minus dataframes
start_index = 1
plus_df_adjusted, start_index = adjust_numbering(plus_df, start_index)
minus_df_adjusted, _ = adjust_numbering(minus_df, start_index)

# Save the adjusted dataframes
minus_output_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_minus_ratio_adjusted_numbering.txt'
plus_output_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_plus_ratio_adjusted_numbering.txt'

minus_df_adjusted.to_csv(minus_output_path, sep="\t", index=False, header=False)
plus_df_adjusted.to_csv(plus_output_path, sep="\t", index=False, header=False)

print("Numbering adjusted and files saved:")
print(f"Plus Strand: {plus_output_path}")
print(f"Minus Strand: {minus_output_path}")
