import numpy as np
import pandas as pd
import os, sys

df = pd.DataFrame(columns = ['accession_name'])
for file in os.listdir('../Transformations/Out'):
    if file[-4:] == '.out':
        tdf = pd.read_csv('../Transformations/Out/%s' %file, sep = ' ', header = None, names = ['accession_name', 'fam', file[:-4]])[['accession_name', file[:-4]]]
        df = pd.merge(df, tdf, on = 'accession_name', how = 'outer')
        
        
#extra
keep = []
for c in df.columns[1:]:
    if c[0] == 'C' and c[1] != 'h':
        keep.append(c)

df = df.loc[:, ['accession_name'] + keep]
    
df.to_csv('../Transformations/68transformed_%s.csv' %sys.argv[1], index = False)
df = pd.read_csv('../Transformations/68transformed_%s.csv' %sys.argv[1], header = None)
np.savetxt('../Transformations/68transformed_%s.txt' %sys.argv[1], df.fillna('NA').values, fmt = '%s', delimiter="\t")