import pandas as pd

# Load the data
minus_data_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_minus_merged.txt'
plus_data_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_plus_merged.txt'

# Read the files into pandas DataFrames
minus_df = pd.read_csv(minus_data_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])
plus_df = pd.read_csv(plus_data_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])

# Combine the 'Value' columns from both DataFrames to compute the total sum
total_sum = minus_df['Value'].sum() + plus_df['Value'].sum()

# Normalize the 'Value' columns by dividing by the total sum
minus_df['Value'] = round((minus_df['Value']*16569) / total_sum, 3)
plus_df['Value'] = round((plus_df['Value']*16569) / total_sum, 3)

# Define a function to adjust the numbering
def adjust_numbering(df):
    # Create contiguous numbering for the 'Start' and 'End' columns
    df['Start'] = range(1, len(df) + 1)
    df['End'] = range(1, len(df) + 1)
    return df

# Apply the function to both DataFrames
minus_df_adjusted = adjust_numbering(minus_df)
plus_df_adjusted = adjust_numbering(plus_df)

# Save the adjusted DataFrames back to text files
minus_output_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_minus_normalized.txt'
plus_output_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_plus_normalized.txt'

minus_df_adjusted.to_csv(minus_output_path, sep="\t", index=False, header=False)
plus_df_adjusted.to_csv(plus_output_path, sep="\t", index=False, header=False)

print("Normalization and adjustment complete. Files saved:")
print(f"Minus Strand: {minus_output_path}")
print(f"Plus Strand: {plus_output_path}")
