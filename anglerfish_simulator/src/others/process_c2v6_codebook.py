import pandas as pd
import numpy as np

c2v6_codebook = pd.read_csv('codebook.tsv', sep='\t')


barcode = np.array(c2v6_codebook.iloc[:,1:])

barcode_str = pd.DataFrame([])

for r in range(barcode.shape[0]):
    
    barcode_str = pd.concat([barcode_str, pd.DataFrame(['  '.join(str(x) for x in barcode[r])])])
    
    
barcode_str.columns = ['barcode']

barcode_str = barcode_str.set_axis([i for i in range(c2v6_codebook.shape[0])], axis='index')

C2V6_codebook_no_distribution = pd.concat([c2v6_codebook.iloc[:,[0]], barcode_str], axis = 1)

C2V6_codebook_no_distribution.to_csv('codebook_harvard_test.csv')