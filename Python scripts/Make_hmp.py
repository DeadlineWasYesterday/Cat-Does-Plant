#converts .vcf to .hmp.txt format for GAPIT
#can be run directly
#example command 'python Make_hmp.py ../Data/WF2.vcf ../Data/WF2'
#note that first arg has the .vcf file extension, second part does not need one.

import numpy as np
import pandas as pd
import sys, io 


def read_vcf(path):
    with open(path, 'r') as f:
        lines = []
        for line in f:
            if line[:2] != '##':
                lines.append(line.replace('|', '/'))


    lines = ''.join(lines)
    return pd.read_csv(io.StringIO(lines), dtype = str, sep='\t')


if __name__ == '__main__':

    df = read_vcf(sys.argv[1])
    
    print('Working with %s' % sys.argv[1])

    t = []
    for i,r in df.iterrows():
        t.append(r.apply(lambda x: x.replace('0', r['REF']).replace('/', '').replace('1', r['ALT'])))

    df3 = pd.DataFrame(t)

    df2 = pd.DataFrame()
    df2['rs'] = df['ID']   
    df2['alleles'] = df['REF'] + '/' + df['ALT']
    df2['chrom'] = df['#CHROM']
    df2['pos'] = df['POS']
    df2['strand'] = np.nan
    df2['assembly'] = np.nan
    df2['center'] = np.nan
    df2['protLSID'] = np.nan
    df2['assayLSID'] = np.nan
    df2['panel'] = np.nan
    df2['QCcode'] = np.nan

    for c in df3.columns[9:]:
        df2[c] = df3[c]

    df2['chrom'] = [int(s[-2:]) for s in df2['chrom'].to_list()]
    l = []
    for a in df2['alleles'].to_list():
        t = a.split('/')
        t.sort()
        l.append('/'.join(t))
    df2['alleles'] = l
    
    print('Columns:')
    print(df2.columns)
    print('Head:')
    print(df2.head())
    print('Shape: %dx%d' %df2.shape)

    df2.to_csv('%s.hmp.csv' %sys.argv[2], sep = '\t', index = False)

    df2 = pd.read_csv('%s.hmp.csv' %sys.argv[2], sep = '\t', header = None, dtype = str)
    np.savetxt('%s.hmp.txt' %sys.argv[2], df2.fillna('NA').values, fmt = '%s', delimiter="\t")
    
    

