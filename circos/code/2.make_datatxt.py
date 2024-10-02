# this don't seem to be working





# Load the necessary libraries
import pandas as pd

# Load the data
minus_data_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_WT/graph/minus_merged.txt'
plus_data_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_plus_ratio_merged.txt'

# Read the files into pandas dataframes, assuming they are space-delimited
minus_df = pd.read_csv(minus_data_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])
plus_df = pd.read_csv(plus_data_path, delim_whitespace=True, header=None, names=["Gene", "Start", "End", "Value"])

# Define a function to adjust the numbering
def adjust_numbering(df):
    # Create a contiguous numbering for the 'Start' and 'End' columns
    df['Start'] = range(1, len(df) + 1)
    df['End'] = range(1, len(df) + 1)
    return df

# Apply the function to both dataframes
minus_df_adjusted = adjust_numbering(minus_df)
plus_df_adjusted = adjust_numbering(plus_df)

# Save the adjusted dataframes back to text files
minus_output_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_WT/graph/minus_adjusted.txt'
plus_output_path = '/Users/ahramahn/Library/CloudStorage/OneDrive-UniversityofMiami/1. Research/DMS-MaP-Seq/FASTKD5/FASTKD_project_KO/graph/KO_plus_ratio_adjusted_numbering.txt'

minus_df_adjusted.to_csv(minus_output_path, sep="\t", index=False, header=False)
plus_df_adjusted.to_csv(plus_output_path, sep="\t", index=False, header=False)
