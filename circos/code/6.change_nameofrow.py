import pandas as pd
import numpy as np

# Paths to the karyotype and data files
karyotype_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/karyotype.human.mt.txt'
plus_data_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_plus_ratio_adjusted_numbering.txt'
minus_data_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_minus_ratio_adjusted_numbering.txt'

# Read the karyotype data
karyotype_cols = ['chr', 'dash', 'name', 'abbrev', 'start', 'end', 'chromosome']
karyotype = pd.read_csv(karyotype_path, sep='\s+', header=None, names=karyotype_cols)

# Function to map a position to the karyotype name
def map_position_to_name(pos, karyotype_df):
    # Find the karyotype row where pos is between start and end
    matched = karyotype_df[(karyotype_df['start'] <= pos) & (pos <= karyotype_df['end'])]
    if not matched.empty:
        return matched.iloc[0]['name']
    else:
        return 'Unknown'

# List of genes to set values to zero
genes_to_zero = ['RNR1', 'RNR2', 'ND6', 'NCR', 'NCR1', 'NCR2', 'un1', 'un2', 'un3', 'un4']

# Identify tRNA genes (any gene name starting with 'tRNA')
tRNA_genes = karyotype[karyotype['name'].str.startswith('tRNA')]['name'].unique()

# Combine all genes to be set to zero
genes_to_zero = genes_to_zero + list(tRNA_genes)

# Process Plus Data
print("Processing Plus Data...")
plus_data = pd.read_csv(plus_data_path, sep='\t', header=None, names=['Gene', 'Start', 'End', 'Value'])

# Map positions to gene names
plus_data['Mapped_Gene'] = plus_data['Start'].apply(lambda x: map_position_to_name(x, karyotype))

# Replace the original gene names with the mapped ones
plus_data['Gene'] = plus_data['Mapped_Gene']

# Drop the temporary 'Mapped_Gene' column
plus_data = plus_data.drop(columns=['Mapped_Gene'])

# Set values to zero for specified genes
plus_data.loc[plus_data['Gene'].isin(genes_to_zero), 'Value'] = 0.0

# Add 0.001 to all values to avoid division by zero in future calculations
plus_data['Value'] += 0.001

# Round the 'Value' column to 3 decimal points
plus_data['Value'] = plus_data['Value'].round(3)

# Save the modified plus data
modified_plus_data_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_plus_data_modified_final.txt'
plus_data.to_csv(modified_plus_data_path, sep='\t', index=False, header=False)
print(f"Modified plus data saved to: {modified_plus_data_path}")

# Process Minus Data
print("Processing Minus Data...")
minus_data = pd.read_csv(minus_data_path, sep='\t', header=None, names=['Gene', 'Start', 'End', 'Value'])

# Map positions to gene names
minus_data['Mapped_Gene'] = minus_data['Start'].apply(lambda x: map_position_to_name(x, karyotype))

# Replace the original gene names with the mapped ones
minus_data['Gene'] = minus_data['Mapped_Gene']

# Drop the temporary 'Mapped_Gene' column
minus_data = minus_data.drop(columns=['Mapped_Gene'])

# Set values to zero for specified genes
minus_data.loc[minus_data['Gene'].isin(genes_to_zero), 'Value'] = 0.0

# Add 0.001 to all values to avoid division by zero in future calculations
minus_data['Value'] += 0.001

# Round the 'Value' column to 3 decimal points
minus_data['Value'] = minus_data['Value'].round(3)

# Save the modified minus data
modified_minus_data_path = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_minus_data_modified_final.txt'
minus_data.to_csv(modified_minus_data_path, sep='\t', index=False, header=False)
print(f"Modified minus data saved to: {modified_minus_data_path}")

print("Data processing complete.")
