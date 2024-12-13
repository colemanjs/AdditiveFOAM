import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_errorbars_from_csv(directory):
    """
    Reads CSV files from a directory and generates an errorbar plot.

    Args:
        directory (str): Path to the directory containing CSV files.

    Returns:
        None
    """
    layer_numbers = []
    means = []
    errors = []

    # Iterate through all files in the directory
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".csv"):
            layer_number = int(''.join(filter(str.isdigit, filename)))
            file_path = os.path.join(directory, filename)

            # Read the CSV file
            data = pd.read_csv(file_path)
            
            # Filter data where 'solid' is 1
            solid_data = data[data['solid'] == 1]
            
            # Compute mean and standard deviation for the 'T' column of filtered data
            mean = solid_data['T'].mean()
            std = solid_data['T'].std()

            # Store values
            layer_numbers.append(layer_number)
            means.append(mean)
            errors.append(std)

    # Plot the errorbar plot
    means = np.array(means)
    errors = np.array(errors)
    plt.xlim(min(layer_numbers), max(layer_numbers))
    
    data_range = max(means + errors) - min(means - errors)
    padding = data_range * 0.1
    plt.ylim(min(means - errors) - padding, max(means + errors) + padding)

    # Plot the errorbar plot
    plt.errorbar(layer_numbers, means, yerr=errors, fmt='o', color='C0', capsize=5)
    plt.plot(layer_numbers, means, linestyle='-', marker='', color='C0')


    plt.xlabel("Layer No.")
    plt.ylabel("Temperature (K)")
    plt.grid(True)
    plt.show()

def main():
    parser = argparse.ArgumentParser(
        description="Generate errorbar plot from CSV files in a directory."
    )
    parser.add_argument(
        "-directory", type=str, required=True,
        help="Directory containing the CSV files."
    )
    args = parser.parse_args()

    # Generate the plot
    plot_errorbars_from_csv(args.directory)

if __name__ == "__main__":
    main()

