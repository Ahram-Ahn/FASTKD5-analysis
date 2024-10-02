import pandas as pd

# Load the data
minus_data_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_minus_ratio_merged_test.txt'
plus_data_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_plus_ratio_merged_test.txt'

# Read the files into pandas dataframes
minus_df = pd.read_csv(minus_data_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])
plus_df = pd.read_csv(plus_data_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])

# Adjust numbering separately for plus and minus dataframes
def adjust_numbering(df):
    # Reset the Start and End to start from 1
    df['Start'] = range(1, len(df) + 1)
    df['End'] = df['Start']
    return df

# Adjust numbering for plus and minus dataframes
plus_df_adjusted = adjust_numbering(plus_df)
minus_df_adjusted = adjust_numbering(minus_df)

# Save the adjusted dataframes
minus_output_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_minus_ratio_adjusted_numbering_test.txt'
plus_output_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_plus_ratio_adjusted_numbering_test.txt'

minus_df_adjusted.to_csv(minus_output_path, sep="\t", index=False, header=False)
plus_df_adjusted.to_csv(plus_output_path, sep="\t", index=False, header=False)

print("Numbering adjusted and files saved:")
print(f"Plus Strand: {plus_output_path}")
print(f"Minus Strand: {minus_output_path}")
