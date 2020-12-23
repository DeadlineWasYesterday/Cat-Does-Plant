#supplement to reverse vcfer for restoring heterozygoes.

import pandas as pd
#df1 = pd.read_csv('hqBSNP2.csv', dtype = str)
df2 = pd.read_csv('hqBSNPlel2.csv', dtype = str)
df3 = pd.read_csv('rvcout.csv', dtype = str)
df4 = pd.read_csv('lqBSNPlel.csv', dtype = str)

df4.iloc[:,1:][df4.iloc[:,1:].notnull()] = './.' 

df = pd.concat([df4[['CP']], df4.iloc[:,1:]], axis = 1)

df = df.sort_values('CP').reset_index(drop = True)
df3 = df3.sort_values('CP').reset_index(drop = True)

df3 = df3.combine_first(df)

df2 = df2.loc[:,df3.columns[1:]].sort_values('CP').reset_index(drop = True)

for c in df2.columns[1:]:
    tdf = df3.loc[:,c][df2.loc[:,c] == '0/1']
    
    for i,r in tdf.iteritems():
        df3.loc[i,c] = '0/' + r.split('/')[1]
    print('%d Done.' %list(df2.columns[1:]).index(c))
    
    
df3.to_csv('rvcwlq.csv', index = False)