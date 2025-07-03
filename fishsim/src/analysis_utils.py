
import gzip
import csv
import pandas as pd
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.stats as stats
from adjustText import adjust_text
import numpy as np


def convert_gz_to_csv(input_filename, output_filename):

    # Open the gzipped TSV file for reading
    with gzip.open(input_filename, 'rt', encoding='utf-8') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        
        # Open the CSV file for writing
        with open(output_filename, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            
            for row in reader:
                writer.writerow(row)
    
    print(f"Converted {input_filename} to {output_filename}")
    
    

def savannah_results_to_bulk(input_filename, output_filename):
    
    ''' Input file is expected to be in csv format.'''
    savannah_spots_results = pd.read_csv(input_filename)
    
    savannah_bulk = savannah_spots_results['target'].value_counts().reset_index()
    
    savannah_bulk.columns = ['gene_name', 'count']
    
    savannah_bulk.to_csv(output_filename, index=False)



def log_correlation(df, plot_name, output_path, fig_name=None):

    f1, ax = plt.subplots(figsize=(9, 9))
    sns.set_palette("deep")
    sns.scatterplot(x="log_tpm",y="log_counts",data=df,ax=ax)
                   
    count_type = "V_counts"
            
        
    ax.set_title("TPM Correlation for " + plot_name, fontsize = 15)
    ax.set_xlabel("log2(TPM+1e-4), V_bulk", fontsize = 20)
    ax.set_ylabel("log2(# detected counts+1), "+ count_type, fontsize = 20)
    
    
    def plotlabel(xvar, yvar, label):
        ax.text(xvar+0.02, yvar, label)
  
        
    pearson, _ = stats.pearsonr(df["log_tpm"],df["log_counts"])
    spearman, _ = stats.spearmanr(df["log_tpm"],df["log_counts"])
    #kendalltau, _ = stats.kendalltau(df["log_tpm"],df["log_counts"])
       
    ax.text(.01, .95, 'Pearson = {:.2f}\nSpearman = {:.2f}'.format(pearson,spearman),transform=ax.transAxes)
    #df.apply(lambda x: plotlabel(x['log_tpm'],  x['log_counts'], x['gene_symbol']), axis=1)
    #plt.axvline(x = 0, color = 'r', label = 'axvline - full height')
    
    texts = []
    for xs,ys,label in zip(df['log_tpm'],df['log_counts'],df['gene_name']):
            texts.append(ax.text(xs,ys,label))
    adjust_text(texts, force_points=0.2, force_text=0.2,expand_points=(1, 1), expand_text=(1, 1), arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.tight_layout()
    
    if fig_name==None:
        plt.savefig(os.path.join(output_path, 'log_correlation.png'))
    else:
        plt.savefig(os.path.join(output_path, f'log_correlation_{fig_name}.png'))

    
    return pearson, spearman




def method_vs_bulk(method_bulk, bulk_seq, plot_name, output_path, fig_name=None):
    
    ''' We assume the method_bulk already has the correct column names: [gene_name, count]'''
    # Remove blanks
    method_bulk = method_bulk[~method_bulk['gene_name'].str.contains('blank_')]    
    method_bulk = method_bulk[~method_bulk['gene_name'].str.contains('Blank-')]    

    
    # Visualize Merlin FOV 18 all Z results
    merged_df = method_bulk.merge(bulk_seq[['gene_symbol', 'bulk_exp']], left_on='gene_name', right_on='gene_symbol', how='left')
    
    # Drop the redundant column
    merged_df = merged_df.drop(columns=['gene_symbol'])
    
    merged_df['log_counts'] = np.log2(merged_df['count']+1)
    
    merged_df['log_tpm'] = np.log2(merged_df['bulk_exp']+0.0001)
    
    pearson_corr, spearman_corr =  log_correlation(merged_df, plot_name, output_path, fig_name)
    
    return pearson_corr, spearman_corr



def method_vs_method(method_1_bulk, method_2_bulk, plot_name, output_path, fig_name=None):
    
    ''' We assume the method_bulk already has the correct column names: [gene_name, count]'''
    # Remove blanks
    method_1_bulk = method_1_bulk[~method_1_bulk['gene_name'].str.contains('blank_')]    
    
    method_2_bulk = method_2_bulk[~method_2_bulk['gene_name'].str.contains('blank_')]    
    
    method_1_bulk = method_1_bulk[~method_1_bulk['gene_name'].str.contains('Blank-')]    
    
    method_2_bulk = method_2_bulk[~method_2_bulk['gene_name'].str.contains('Blank-')] 

    method_1_bulk = method_1_bulk.rename(columns={'count': 'method_1_count'})

    method_2_bulk = method_2_bulk.rename(columns={'count': 'method_2_count'})

    # Visualize Merlin FOV 18 all Z results
    merged_df = method_1_bulk.merge(method_2_bulk[['gene_name', 'method_2_count']], left_on='gene_name', right_on='gene_name', how='left')
    
    # Drop the redundant column    
    merged_df['log_counts_method_1'] = np.log2(merged_df['method_1_count']+1)
    
    merged_df['log_counts_method_2'] = np.log2(merged_df['method_2_count']+1)
    
    pearson_corr, spearman_corr =  log_correlation_mvm(merged_df, plot_name, output_path, fig_name)
    
    return pearson_corr, spearman_corr



def log_correlation_mvm(df, plot_name, output_path, fig_name=None):

    f1, ax = plt.subplots(figsize=(9, 9))
    sns.set_palette("deep")
    sns.scatterplot(x="log_counts_method_1",y="log_counts_method_2",data=df,ax=ax)
                   
    count_type = "V_counts"
            
        
    ax.set_title("Log Correlation for " + plot_name, fontsize = 15)
    ax.set_xlabel("log2(# detected counts+1), " + count_type, fontsize = 20)
    ax.set_ylabel("log2(# detected counts+1), " + count_type, fontsize = 20)
    
    
    def plotlabel(xvar, yvar, label):
        ax.text(xvar+0.02, yvar, label)
  
        
    pearson, _ = stats.pearsonr(df["log_counts_method_1"],df["log_counts_method_2"])
    spearman, _ = stats.spearmanr(df["log_counts_method_1"],df["log_counts_method_2"])
    #kendalltau, _ = stats.kendalltau(df["log_tpm"],df["log_counts"])
       
    ax.text(.01, .95, 'Pearson = {:.2f}\nSpearman = {:.2f}'.format(pearson,spearman),transform=ax.transAxes)
    #df.apply(lambda x: plotlabel(x['log_tpm'],  x['log_counts'], x['gene_symbol']), axis=1)
    #plt.axvline(x = 0, color = 'r', label = 'axvline - full height')
    
    texts = []
    for xs,ys,label in zip(df['log_counts_method_1'],df['log_counts_method_2'],df['gene_name']):
            texts.append(ax.text(xs,ys,label))
    adjust_text(texts, force_points=0.2, force_text=0.2,expand_points=(1, 1), expand_text=(1, 1), arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.tight_layout()
    
    if fig_name==None:
        plt.savefig(os.path.join(output_path, 'log_correlation.png'))
    else:
        plt.savefig(os.path.join(output_path, f'log_correlation_{fig_name}.png'))
    
    return pearson, spearman




    
def merlin_sv_to_bulk(input_merlin_sv, codebook_tsv_path, output_filepath, fov=None):

    if fov != None:
        
        input_merlin_sv = input_merlin_sv[(input_merlin_sv['fov'] == fov)]

    input_merlin_sv = input_merlin_sv['barcode_id'].value_counts().reset_index()

    input_merlin_sv.columns = ['barcode_id', 'count']

    input_merlin_sv = input_merlin_sv.sort_values(by='barcode_id').reset_index(drop=True)

    # Load the codebook's target for matching
    codebook = pd.read_csv(codebook_tsv_path, sep='\t').reset_index()[['index', 'target']]

    input_merlin_sv_bulk = input_merlin_sv.merge(codebook, left_on='barcode_id', right_on='index', how='left')

    input_merlin_sv_bulk.drop('index', axis=1, inplace=True)
    
    input_merlin_sv_bulk = input_merlin_sv_bulk[['barcode_id', 'target', 'count']]
    
    input_merlin_sv_bulk.columns = ['barcode_id', 'gene_name', 'count']

    # Remove the blanks
    input_merlin_sv_bulk = input_merlin_sv_bulk[~input_merlin_sv_bulk['gene_name'].str.contains('blank_')]    

    input_merlin_sv_bulk.to_csv(output_filepath)
    
    return input_merlin_sv_bulk
    