import pandas as pd
import numpy as np
import os

# Define file paths for WT and KO samples
file_paths = {
    'WT': '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/WT_plus_data_modified_final.txt',
    'KO': '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/KO_plus_data_modified_final.txt',
}

# Define output directory for Circos data files
output_dir = '/Users/ahramahn/Documents/circos-0.69-9/3.delta_DMS/WT_vs_KO_output_circos_data/'
os.makedirs(output_dir, exist_ok=True)

# Parameters
window_size = 5  # Window size for smoothing

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
    Ensures that all positions are included for each unit.
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
        min_start = int(min(unit_wt['start'].min() if not unit_wt.empty else np.inf,
                            unit_ko['start'].min() if not unit_ko.empty else np.inf))
        max_end = int(max(unit_wt['end'].max() if not unit_wt.empty else -np.inf,
                          unit_ko['end'].max() if not unit_ko.empty else -np.inf))
        
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

def get_combined_top95_percentile(wt_values, ko_values):
    """
    Calculates the 95th percentile from the combined non-zero values of WT and KO.
    This is used for normalization.
    """
    combined_non_zero = pd.concat([wt_values, ko_values])
    non_zero_values = combined_non_zero[combined_non_zero > 0]
    if len(non_zero_values) == 0:
        return 1  # Avoid division by zero if all values are zero
    return np.nanpercentile(non_zero_values, 95)

def scale_values(values, scale_factor):
    """
    Scales the values by the given scaling factor.
    Replaces any zero values with a small number to avoid division by zero.
    """
    # Replace zero with a small number (e.g., 0.001) to avoid division by zero
    scaled = values.replace(0, 0.001) / scale_factor
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
    """
    # Sort data by unit and start position
    data = data.sort_values(['unit', 'start']).reset_index(drop=True)
    
    # Initialize a list to collect binned data
    binned_list = []
    
    for unit in data['unit'].unique():
        unit_data = data[data['unit'] == unit]
        num_bins = int(np.ceil(len(unit_data) / window_size))
        
        for i in range(num_bins):
            bin_start_idx = i * window_size
            bin_end_idx = min((i + 1) * window_size, len(unit_data))
            bin_data = unit_data.iloc[bin_start_idx:bin_end_idx]
            
            bin_start = bin_data['start'].min()
            bin_end = bin_data['end'].max()
            bin_mean_log2_ratio = bin_data['log2_ratio'].mean()
            
            binned_list.append({
                'unit': unit,
                'start': bin_start,
                'end': bin_end,
                'log2_ratio': bin_mean_log2_ratio
            })
    
    binned_df = pd.DataFrame(binned_list)
    return binned_df

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
    
    # Ensure values are numeric
    wt_data['value'] = pd.to_numeric(wt_data['value'], errors='coerce').fillna(0)
    ko_data['value'] = pd.to_numeric(ko_data['value'], errors='coerce').fillna(0)
    
    # Get complete positions
    complete_df = get_complete_positions(wt_data, ko_data)
    
    # Compute combined scaling factor (95th percentile) from both WT and KO
    combined_top95_percentile = get_combined_top95_percentile(complete_df['value_WT'], complete_df['value_KO'])
    print(f"Combined 95th percentile for scaling: {combined_top95_percentile}")
    
    # Scale WT and KO values using the combined scaling factor
    complete_df['scaled_WT'] = scale_values(complete_df['value_WT'], combined_top95_percentile)
    complete_df['scaled_KO'] = scale_values(complete_df['value_KO'], combined_top95_percentile)
    
    # Compute log2 ratio
    complete_df['log2_ratio'] = compute_log2_ratio(complete_df['scaled_KO'], complete_df['scaled_WT'])
    
    # Handle infinite or undefined log2 ratios
    complete_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    complete_df['log2_ratio'] = complete_df['log2_ratio'].fillna(0)
    
    # Round log2_ratio to three decimal points
    complete_df['log2_ratio'] = complete_df['log2_ratio'].round(3)
    
    # Save single nucleotide resolution log2 ratios
    snr_output_file = os.path.join(output_dir, "WT_vs_KO_log2_ratio_single_nucleotide.txt")
    complete_df[['unit', 'start', 'end', 'log2_ratio']].to_csv(
        snr_output_file, sep='\t', header=False, index=False, float_format='%.3f')
    print(f"Single nucleotide log2 ratios saved to {snr_output_file}")
    
    # Apply binning/smoothing
    binned_df = bin_data(complete_df, window_size)
    
    # Round log2_ratio to three decimal points in binned data
    binned_df['log2_ratio'] = binned_df['log2_ratio'].round(3)
    
    # Save binned log2 ratios
    binned_output_file = os.path.join(output_dir, f"WT_vs_KO_log2_ratio_binned_window_{window_size}.txt")
    binned_df[['unit', 'start', 'end', 'log2_ratio']].to_csv(
        binned_output_file, sep='\t', header=False, index=False, float_format='%.3f')
    print(f"Binned log2 ratios saved to {binned_output_file}\n")
    
    print("Processing complete.")

# Execute the processing
process_samples(file_paths, window_size, output_dir)
