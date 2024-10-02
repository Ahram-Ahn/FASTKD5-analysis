import pandas as pd
import numpy as np

# Paths to the modified KO and WT plus data
KO_plus_data_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_plus_data_modified_final.txt'
WT_plus_data_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/WT_plus_data_modified_final.txt'  # Ensure you have this file

# Read the data
KO_plus_df = pd.read_csv(KO_plus_data_path, sep='\t', header=None, names=['Gene', 'Start', 'End', 'KO_Value'])
WT_plus_df = pd.read_csv(WT_plus_data_path, sep='\t', header=None, names=['Gene', 'Start', 'End', 'WT_Value'])

# Merge the data on 'Gene', 'Start', and 'End'
merged_plus_df = pd.merge(KO_plus_df, WT_plus_df, on=['Gene', 'Start', 'End'])

# Compute log2 fold difference
merged_plus_df['Log2_Fold_Change'] = np.log2(merged_plus_df['KO_Value'] / merged_plus_df['WT_Value'])

# Round the log fold change to three decimal points
merged_plus_df['Log2_Fold_Change'] = merged_plus_df['Log2_Fold_Change'].round(3)

# Save the final data
final_plus_output_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/WT_vs_KO_log_fold_difference_plu_final.txt'
merged_plus_df[['Gene', 'Start', 'End', 'Log2_Fold_Change']].to_csv(final_plus_output_path, sep='\t', index=False, header=False)

print(f"Log fold difference for plus data computed and saved to: {final_plus_output_path}")
