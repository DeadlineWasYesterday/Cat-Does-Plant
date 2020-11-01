import numpy as np
import pandas as pd
import io

df = pd.read_csv('WF2.vcf', dtype = str, sep = '\t')

dfind = df.loc[:, ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
          'IRIS_313-10340',
 'IRIS_313-10341',
 'IRIS_313-10539',
 'IRIS_313-10971',
 'IRIS_313-10972',
 'IRIS_313-10973',
 'IRIS_313-10974',
 'IRIS_313-10975',
 'IRIS_313-10977',
 'IRIS_313-10978',
 'IRIS_313-10982',
 'IRIS_313-10983',
 'IRIS_313-10985',
 'IRIS_313-10986',
 'IRIS_313-11113',
 'IRIS_313-11114',
 'IRIS_313-11115',
 'IRIS_313-11203',
 'IRIS_313-11204',
 'IRIS_313-11205',
 'IRIS_313-11206',
 'IRIS_313-11207',
 'IRIS_313-11208',
 'IRIS_313-11211',
 'IRIS_313-11212',
 'IRIS_313-11219',
 'IRIS_313-11220',
 'IRIS_313-11221',
 'IRIS_313-11222',
 'IRIS_313-11224',
 'IRIS_313-11225',
 'IRIS_313-11226',
 'IRIS_313-11227',
 'IRIS_313-11228',
 'IRIS_313-11229',
 'IRIS_313-11230',
 'IRIS_313-11231',
 'IRIS_313-11241',
 'IRIS_313-11321',
 'IRIS_313-11326',
 'IRIS_313-11399',
 'IRIS_313-11400',
 'IRIS_313-11402',
 'IRIS_313-11404',
 'IRIS_313-11487',
 'IRIS_313-11558',
 'IRIS_313-11722',
 'IRIS_313-11945',
 'IRIS_313-11964',
 'IRIS_313-8244',
 'IRIS_313-8349',
 'IRIS_313-8437',
 'IRIS_313-8509',
 'IRIS_313-8632',
 'IRIS_313-8683',
 'IRIS_313-8703',
 'IRIS_313-8717',
 'IRIS_313-8733',
 'IRIS_313-8737',
 'IRIS_313-8850',
 'IRIS_313-8854',
 'IRIS_313-8930',
 'IRIS_313-8932',
 'IRIS_313-9049',
 'IRIS_313-9066',
 'IRIS_313-9067',
 'IRIS_313-9072',
 'IRIS_313-9108',
 'IRIS_313-9139',
 'IRIS_313-9148',
 'IRIS_313-9156',
 'IRIS_313-9174',
 'IRIS_313-9218',
 'IRIS_313-9262',
 'IRIS_313-9271',
 'IRIS_313-9325',
 'IRIS_313-9391',
 'IRIS_313-9397',
 'IRIS_313-9594',
 'IRIS_313-9606',
 'IRIS_313-9617',
 'IRIS_313-9696']]

#ref allele filter
keep = []
for i,r in dfind.iterrows():
    if ''.join(r).count('0/0') > 4:
        keep.append(i)
        
dfind = dfind.loc[keep, :].reset_index(drop = True)
dfind.to_csv('WF2indR4.vcf', index = False, sep = '\t')

keep = []
for i,r in dfind.iterrows():
    if ''.join(r).count('0/0') > 7:
        keep.append(i)
        
dfind = dfind.loc[keep, :]
dfind.to_csv('WF2indR7.vcf', index = False, sep = '\t')


