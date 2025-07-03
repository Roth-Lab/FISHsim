"""
Given bulk-seq with TPM/FPKM and gene symbol and the target experiment codebook, we generate a 'special codebook' with the dedicated simulated counts according to the bulk-seq.

This special codebook is used to simulate testing data to test the trained anglerFISH model.
"""

import pandas as pd
import numpy as np
import os



def generate_experiment_codebook(bulk_seq_path, codebook_path, output_path, target_total_emitter_count):
    
    """ 
    codebook path: csv format
    bulk seq path: csv.gz format
    target_total_emitter_count: should match the emitter count from the simulated data used for training the anglerFISH model
    """
    bulk_seq = pd.read_csv(bulk_seq_path, compression='gzip')

    codebook_gene_names = pd.read_csv(codebook_path, skiprows=3)['name']

    # Find common genes in bulk-seq and codebook
    common_genes = bulk_seq[bulk_seq["gene_symbol"].isin(codebook_gene_names)]

    # Return the relevant columns
    result = common_genes[["gene_symbol", "bulk_exp"]]

    result['bulk_exp_proportion'] = result['bulk_exp']/result['bulk_exp'].sum()

    result['simulated_count'] = np.floor(target_total_emitter_count*result['bulk_exp_proportion']).astype(np.int16)
    
    # Set the index of the result DataFrame to 'gene_symbol'
    result.set_index('gene_symbol', inplace=True)
    
    # Reorder the result DataFrame based on the order of codebook['target']
    reordered_result = result.reindex(codebook_gene_names).reset_index()

    new_codebook = pd.read_csv(codebook_path)
    
    new_codebook.loc[2, 'Unnamed: 4'] = 'distribution'
    
    new_codebook.loc[3:, 'Unnamed: 4'] = reordered_result['simulated_count'].values
    
    new_codebook.to_csv(output_path)



if __name__ == "__main__":

    
    data_directory = 'C:/Users/jenkints/Documents/GitHub/savannah_categorical/data/xp8054'
    
    anglerfish_simulator_directory = 'C:/Users/jenkints/Documents/GitHub/anglerfish_simulator/anglerfish_simulator'
    
    bulk_seq_path = os.path.join(data_directory, 'XP2474_4T1_C1E1_new_bulk.csv.gz')

    codebook_path = os.path.join(anglerfish_simulator_directory, 'resources', 'codebook_xp7049.csv')
    
    output_path = os.path.join(anglerfish_simulator_directory, 'resources', 'codebook_xp8054_distribution.csv')
    
    target_total_emitter_count = 40000
    
    generate_experiment_codebook(bulk_seq_path, codebook_path, output_path, target_total_emitter_count)