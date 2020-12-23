#After warpedlmm has done its thing, run to compile the results into one file.

import numpy as np
import pandas as pd
import os, sys

df = pd.DataFrame(columns = ['accession_name'])
for file in os.listdir('../Transformations/Out'):
    if file[-4:] == '.out':
        tdf = pd.read_csv('../Transformations/Out/%s' %file, sep = ' ', header = None, names = ['accession_name', 'fam', file[:-4]])[['accession_name', file[:-4]]]
        df = pd.merge(df, tdf, on = 'accession_name', how = 'outer')
    
df.to_csv('../Transformations/39transformed_%s.csv' %sys.argv[1], index = False)
df = pd.read_csv('../Transformations/39transformed_%s.csv' %sys.argv[1], header = None)
np.savetxt('../Transformations/39transformed_%s.txt' %sys.argv[1], df.fillna('NA').values, fmt = '%s', delimiter="\t")