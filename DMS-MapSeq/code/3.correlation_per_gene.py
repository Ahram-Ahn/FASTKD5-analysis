import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import r2_score
import os

# Load WT and KO data
wt_file = '/Users/ahramahn/Documents/circos-0.69-9/code_for_delta_dms/WT_plus_data_modified_final.txt'
ko_file = '/Users/ahramahn/Documents/circos-0.69-9/code_for_delta_dms/KO_plus_data_modified_final.txt'

# Column names based on provided data format
columns = ['gene', 'position', 'index', 'mutation_ratio']

# Read the files
wt_data = pd.read_csv(wt_file, sep='\t', names=columns)
ko_data = pd.read_csv(ko_file, sep='\t', names=columns)

# Replace "ATP8/6" with "ATP86" to avoid errors
wt_data['gene'] = wt_data['gene'].replace('ATP8/6', 'ATP86')
ko_data['gene'] = ko_data['gene'].replace('ATP8/6', 'ATP86')

# Merge WT and KO data on gene and position to compare
merged_data = pd.merge(wt_data, ko_data, on=['gene', 'position'], suffixes=('_WT', '_KO'))

# Remove rows where either WT or KO has NaN values
merged_data = merged_data.dropna(subset=['mutation_ratio_WT', 'mutation_ratio_KO'])

# Filter out rows where both WT and KO mutation ratios are 0.001 and higher than 0.15
filtered_data = merged_data[
    ((merged_data['mutation_ratio_WT'] > 0.001) & (merged_data['mutation_ratio_KO'] > 0.001)) &
    (merged_data['mutation_ratio_WT'] <= 0.15) & (merged_data['mutation_ratio_KO'] <= 0.15)
]

# Create a directory to save the plots
os.makedirs('correlation_plots', exist_ok=True)

# Prepare list to store correlation data for CSV
correlation_data = []

# Set the font to Arial
plt.rcParams["font.family"] = "Arial"

# Iterate over each gene and plot
unique_genes = filtered_data['gene'].unique()
for gene in unique_genes:
    gene_data = filtered_data[filtered_data['gene'] == gene]
    
    # Scatter plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='mutation_ratio_WT', y='mutation_ratio_KO', data=gene_data)

    # Remove X and Y axis labels and ticks
    plt.xlabel('')
    plt.ylabel('')


    # Calculate R-squared
    if len(gene_data) > 1:  # Ensure there's enough data to calculate R-squared
        r2 = r2_score(gene_data['mutation_ratio_WT'], gene_data['mutation_ratio_KO'])
    else:
        r2 = None  # If not enough data points, set R-squared to None

    # Save correlation data for CSV file
    correlation_data.append({
        'gene': gene,
        'num_positions': len(gene_data),
        'r_squared': r2
    })
    
    # Annotate plot with R-squared value, handling None case
    if r2 is not None:
        plt.title(f'{gene} (R² = {r2:.4f})', fontweight='bold', fontsize=14)
    else:
        plt.title(f'{gene} (R² = N/A)', fontweight='bold', fontsize=14)
    
    # Clean up gene name for saving the plot (remove any slashes or problematic characters)
    safe_gene_name = gene.replace('/', '_')
    
    # Save the plot as both PNG and SVG formats
    plt.savefig(f'correlation_plots/{safe_gene_name}_correlation.png')
    plt.savefig(f'correlation_plots/{safe_gene_name}_correlation.svg', format='svg')
    plt.close()

# Save the correlation data to a CSV file
correlation_df = pd.DataFrame(correlation_data)
correlation_df.to_csv('plus_correlation_data.csv', index=False)

print("Correlation plots saved in 'plus_correlation_plots' directory as PNG and SVG.")
print("Correlation data saved to 'plus_correlation_data.csv'.")
