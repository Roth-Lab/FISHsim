import pandas as pd
import numpy as np

def generate_gene_distributions(input_path="../codebooks/C1E1_codebook_no_distribution.csv",
                                targeted_emitter_count=5000):
    """
    Generates gene distributions based on Dirichlet and predefined distributions,
    then scales them so that the total sum of the "distribution" column equals targeted_emitter_count.
    Saves the outputs as CSV files.

    Parameters:
    - input_path (str): Path to the input CSV file.
    - targeted_emitter_count (int): Total sum of the "distribution" column in the output.

    Output:
    - Saves three CSV files with different distributions.
    """

    # Load the full dataset
    df_full = pd.read_csv(input_path, header=None)

    # Extract headers (first three rows)
    headers = df_full.iloc[:3]

    # Extract main data (starting from row 4)
    df = df_full.iloc[3:].reset_index(drop=True)

    # Set the proper column names from the fourth row
    df.columns = df.iloc[0]  # Fourth row becomes column names
    df = df[1:].reset_index(drop=True)  # Remove the row used for column names

    # Keep only necessary columns
    df = df[['name', 'numeric_id', 'id', 'barcode']]

    # Define distribution types
    num_genes = len(df)
    distribution_types = {
        "random_distribution": np.random.dirichlet(np.ones(num_genes) * 0.5),  # Random varying distribution
        "uniform_distribution": np.random.dirichlet(np.ones(num_genes) * 10),
        "extreme_distribution": np.random.dirichlet(np.ones(num_genes) * 0.1)
    }

    # Assign extreme distribution to Cd4
    if "Cd4" in df["name"].values:
        idx_cd4 = df[df["name"] == "Cd4"].index[0]
        distribution_types["extreme_distribution"][idx_cd4] = 1000  # Extreme value for Cd4

    # Scale distributions so they sum up to targeted_emitter_count
    for key in distribution_types:
        distribution_types[key] *= targeted_emitter_count / np.sum(distribution_types[key])
        distribution_types[key] = np.round(distribution_types[key]).astype(int)
        
        # Adjust rounding error to maintain exact sum
        # This only applies to extreme_distribution and random_distribution, as issues will arise
        # with uniform_distribution due to having a non-zero remainder issue
        if key != "uniform_distribution":
            diff = targeted_emitter_count - np.sum(distribution_types[key])
            if diff != 0:
                max_idx = np.argmax(distribution_types[key])  # Adjust the highest value to fix sum
                distribution_types[key][max_idx] += diff


    # Generate and save CSVs while preserving headers
    output_files = []
    for dist_type, values in distribution_types.items():
        df["distribution"] = values  # Add distribution column

        # Add a new column to the header row for "distribution"
        headers_with_distribution = headers.copy()
        
        # Add 13 columns filled with NaN values
        num_extra_columns = headers_with_distribution.shape[1] - df.shape[1]
        nan_columns = pd.DataFrame(np.nan, index=df.index, columns=[np.nan] * num_extra_columns)
    
        # Concatenate the original dataframe with the NaN columns
        df_expanded = pd.concat([df, nan_columns], axis=1)
        
        # Save the column names
        column_names_df = pd.DataFrame([df_expanded.columns.tolist()])
        
        # Reset column names to unique integer values first to avoid index errors
        df_expanded.columns = range(df_expanded.shape[1])
        
        # Move column names into the first row of the dataframe
        df_expanded = pd.concat([column_names_df, df_expanded], ignore_index=True)

        # Concatenate headers and data
        df_final = pd.concat([headers_with_distribution, df_expanded], ignore_index=True)

        # Save to file
        output_path = f"../codebooks/C1E1_codebook_{dist_type}_emittercount{targeted_emitter_count}.csv"
        df_final.to_csv(output_path, index=False, header=False)
        output_files.append(output_path)

    return output_files



if __name__ == "__main__":

    generate_gene_distributions()
    # Example usage (uncomment to run locally):
    # output_files = generate_gene_distributions("C1E1_codebook_no_distribution.csv")
    # print(output_files)

