import pandas as pd
import numpy as np
import os

# Define file paths for WT and KO samples
file_paths = {
    'WT': '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/plus_ratio_adjusted_numbering_labeled.txt',  # WT data file path
    'KO': '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_plus_ratio_adjusted_numbering_labeled.txt',  # KO data file path
}

# Define output directory for Circos data files
output_dir = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/WT_vs_KO_output_circos_data/'  # Output directory
os.makedirs(output_dir, exist_ok=True)

# Parameters
window_size = 30  # Window size for smoothing

def read_data(file_path):
    """
    Reads the DMS reactivity data from a text file.
    
    Assumes the file has four columns:
    unit, start, end, value
    Separated by tabs and without a header.
    """
    try:
        data = pd.read_csv(file_path, sep='\t', header=None, names=['unit', 'start', 'end', 'value'])
        return data
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def get_complete_positions(wt_data, ko_data):
    """
    Ensures that all positions from 1 to 16569 are included for each unit.
    If a position is missing in either WT or KO data, it is filled with 0.
    """
    # Determine all unique units
    units = pd.concat([wt_data['unit'], ko_data['unit']]).unique()
    
    # Initialize a list to collect complete data
    complete_data = []
    
    for unit in units:
        # Extract data for the current unit
        unit_wt = wt_data[wt_data['unit'] == unit]
        unit_ko = ko_data[ko_data['unit'] == unit]
        
        if unit_wt.empty and unit_ko.empty:
            continue  # Skip if no data for this unit
        
        # Determine the range of positions for the unit
        min_start = min(unit_wt['start'].min() if not unit_wt.empty else np.inf,
                       unit_ko['start'].min() if not unit_ko.empty else np.inf)
        max_end = max(unit_wt['end'].max() if not unit_wt.empty else -np.inf,
                     unit_ko['end'].max() if not unit_ko.empty else -np.inf)
        
        # Handle cases where unit has no WT or KO data
        if min_start == np.inf:
            min_start = unit_ko['start'].min()
        if max_end == -np.inf:
            max_end = unit_ko['end'].max()
        
        # Create a DataFrame with all positions for the unit
        positions = pd.DataFrame({
            'unit': unit,
            'start': range(min_start, max_end + 1),
            'end': range(min_start, max_end + 1)
        })
        
        # Merge with WT and KO data to get the 'value' columns
        positions = positions.merge(unit_wt[['start', 'value']], on='start', how='left').rename(columns={'value': 'value_WT'})
        positions = positions.merge(unit_ko[['start', 'value']], on='start', how='left').rename(columns={'value': 'value_KO'})
        
        # Fill NaNs with 0
        positions['value_WT'] = positions['value_WT'].fillna(0)
        positions['value_KO'] = positions['value_KO'].fillna(0)
        
        complete_data.append(positions)
    
    # Concatenate all units
    complete_df = pd.concat(complete_data, ignore_index=True)
    return complete_df

def get_combined_top10_ratio(wt_values, ko_values):
    """
    Calculates the 95th percentile (top 5%) ratio from the combined non-zero mutation values of WT and KO.
    This is used for normalization.
    """
    combined_non_zero = pd.concat([wt_values, ko_values])
    non_zero_values = combined_non_zero[combined_non_zero > 0]
    if len(non_zero_values) == 0:
        return 1  # Avoid division by zero if all values are zero
    return np.nanpercentile(non_zero_values, 95)

def scale_values(values, scale_ratio):
    """
    Scales the mutation values by the given scaling ratio.
    Replaces any zero values with a small number to avoid division by zero.
    """
    # Replace zero with a small number (e.g., 0.001) to avoid division by zero
    scaled = values.replace(0, 0.001) / scale_ratio
    return scaled

def compute_log2_ratio(ko_scaled, wt_scaled):
    """
    Computes the log2 ratio between KO and WT scaled values.
    """
    return np.log2(ko_scaled / wt_scaled)

def bin_data(data, window_size):
    """
    Applies binning (grouping) to the data based on the window size.
    Groups data into bins of 'window_size' and computes the mean log2 ratio for each bin.
    
    Windows are sequentially grouped as 1-10, 11-20, 21-30, etc., within each unit.
    """
    # Sort data by unit and start position
    data = data.sort_values(['unit', 'start']).reset_index(drop=True)
    
    # Calculate relative positions within each unit
    data['relative_pos'] = data.groupby('unit')['start'].transform(lambda x: x - x.min() + 1)
    
    # Assign bin indices based on relative positions
    data['bin'] = ((data['relative_pos'] - 1) // window_size).astype(int)
    
    # Group by unit and bin, then compute the mean log2 ratio
    binned = data.groupby(['unit', 'bin'])['log2_ratio'].mean().reset_index()
    
    # Calculate bin start and end positions based on relative positions
    bin_positions = data.groupby(['unit', 'bin']).agg({
        'start': 'min',
        'end': 'max'
    }).reset_index()
    
    # Merge bin ratios with bin positions
    binned = binned.merge(bin_positions, on=['unit', 'bin'])
    
    # Optional: Sort the binned data by unit and bin for better readability
    binned = binned.sort_values(['unit', 'bin']).reset_index(drop=True)
    
    return binned[['unit', 'start', 'end', 'log2_ratio']]

def process_samples(file_paths, window_size, output_dir):
    """
    Processes WT and KO samples to compute log2 ratios and outputs Circos-compatible data files.
    """
    # Read WT and KO data
    wt_data = read_data(file_paths['WT'])
    ko_data = read_data(file_paths['KO'])
    
    if wt_data is None or ko_data is None:
        print("WT or KO data could not be read. Exiting.")
        return
    
    # Get complete positions
    complete_df = get_complete_positions(wt_data, ko_data)
    
    # Compute combined scaling ratio (95th percentile) from both WT and KO
    combined_top10_ratio = get_combined_top10_ratio(complete_df['value_WT'], complete_df['value_KO'])
    print(f"Combined top 10% ratio: {combined_top10_ratio}")
    
    # Scale WT and KO values using the combined scaling ratio
    complete_df['scaled_WT'] = scale_values(complete_df['value_WT'], combined_top10_ratio)
    complete_df['scaled_KO'] = scale_values(complete_df['value_KO'], combined_top10_ratio)
    
    # Compute log2 ratio
    complete_df['log2_ratio'] = compute_log2_ratio(complete_df['scaled_KO'], complete_df['scaled_WT'])
    
    # Handle infinite or undefined log2 ratios
    complete_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    complete_df['log2_ratio'] = complete_df['log2_ratio'].fillna(0)  # Define log2(1) = 0 for positions where both KO and WT are zero
    
    # Save single nucleotide resolution log2 ratios (optional)
    snr_output_file = os.path.join(output_dir, "WT_vs_KO_log2_ratio_single_nucleotide.txt")
    complete_df[['unit', 'start', 'end', 'log2_ratio']].to_csv(snr_output_file, sep='\t', header=False, index=False, float_format='%.4f')
    print(f"Single nucleotide log2 ratios saved to {snr_output_file}")
    
    # Apply binning/smoothing
    binned = bin_data(complete_df, window_size)
    
    # Save binned log2 ratios
    binned_output_file = os.path.join(output_dir, f"WT_vs_KO_log2_ratio_binned_window_{window_size}.txt")
    binned.to_csv(binned_output_file, sep='\t', header=False, index=False, float_format='%.4f')
    print(f"Binned log2 ratios saved to {binned_output_file}\n")
    
    # Optionally, you can generate additional plots or perform further analysis here

def plot_log2_ratios(binned_file, window_size, output_plot):
    """
    Plots the binned log2 ratios for visualization.
    """
    import matplotlib.pyplot as plt
    
    # Read the binned data
    binned = pd.read_csv(binned_file, sep='\t', header=None, names=['unit', 'start', 'end', 'log2_ratio'])
    
    # Plot settings
    plt.figure(figsize=(20, 6))
    plt.bar(binned['start'], binned['log2_ratio'], width=window_size, color=np.where(binned['log2_ratio'] >= 0, 'green', 'red'))
    plt.xlabel('Nucleotide Position')
    plt.ylabel('Log2 Ratio (KO/WT)')
    plt.title(f'Binned Log2 Ratio (Window Size = {window_size})')
    plt.axhline(0, color='orange', linestyle='--', linewidth=1.5)
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.show()
    print(f"Plot saved to {output_plot}")

# Execute the processing
process_samples(file_paths, window_size, output_dir)

# Example of plotting the binned data (uncomment if needed)
# binned_file = os.path.join(output_dir, f"WT_vs_KO_log2_ratio_binned_window_{window_size}.txt")
# output_plot = os.path.join(output_dir, f"WT_vs_KO_log2_ratio_binned_window_{window_size}.png")
# plot_log2_ratios(binned_file, window_size, output_plot)

print("All samples processed successfully.")
