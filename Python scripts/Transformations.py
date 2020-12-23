#This script houses the python code for running transformations by WarpedLMM.
#Run transformations_p1 to create the bash script, run it and 
#then run p2 to compile the results.

import numpy as np
import pandas as pd
import os, sys


df = pd.read_csv('../Data/39traits.csv')

print('Read phenotype file.')

for c in df.columns[1:]:
    tdf = df.loc[:, ['accession_name', 'accession_name', c]].dropna()
    tdf.to_csv('../Transformations/In/%s.in' %c, header = None, sep = '\t', index = False)

print('Made input files')
    
for file in os.listdir('../Transformations/In'):
    with open('../Transformations/transform.sh', 'w') as f:
        f.write('#!/bin/bash')
        f.write('\n')
        f.write('conda init bash')
        f.write('\n')
        f.write('conda activate py2')
        f.write('\n')
        f.write('python -m warpedlmm ../Data/%s In/%s --save --output_directory Out' %(sys.argv[1], file))
    os.system('chmod +x ../Transformations/transform.sh')
    os.system('../Transformations/transform.sh')
    os.rename('../Transformations/Out/warpedlmm_pheno.txt', '../Transformations/Out/%s.out' %file[:-3])
    os.rename('../Transformations/Out/warpedlmm_results.txt', '../Transformations/Out/%s_results.out' %file[:-3])
        
print('Ran shell script.')
        
df = pd.DataFrame(columns = ['accession_name'])
for file in os.listdir('../Transformations/Out'):
    tdf = pd.read_csv('../Transformations/%s' %file, header = None, names = ['accession_name', 'fam', file[:-4]])[['accession_name', file[:-4]]]
    df = pd.merge(df, tdf, on = 'accession_name', how = 'outer')
    
df.to_csv('../Transformations/39transformed.csv')
np.savetxt('../Transformations/39transformed.txt', df.fillna('NA').values, fmt = '%s', delimiter="\t")




