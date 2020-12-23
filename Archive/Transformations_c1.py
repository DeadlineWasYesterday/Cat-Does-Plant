#Copy of transformations_p1 for a different file

import numpy as np
import pandas as pd
import os, sys


df = pd.read_csv('../Data/68composite.csv')

print('Read phenotype file.')

for c in df.columns[1:]:
    tdf = df.loc[:, ['accession_name', 'accession_name', c]].dropna()
    tdf.to_csv('../Transformations/In/%s.in' %c, header = None, sep = '\t', index = False)

print('Made input files')

with open('../Transformations/transform.sh', 'w') as f:
    f.write('#!/bin/bash')
    f.write('\n')
    
for file in os.listdir('../Transformations/In'):
    with open('../Transformations/transform.sh', 'a') as f:
        f.write('python -m warpedlmm ../Data/%s In/%s --save --output_directory Out' %(sys.argv[1], file))
        f.write('\n')
        f.write('mv Out/warpedlmm_pheno.txt Out/%s.out' %file[:-3])
        f.write('\n')
        f.write('mv Out/warpedlmm_results.txt Out/%s.results' %file[:-3])
        f.write('\n')
        
        
print('Made shell script.')
        



