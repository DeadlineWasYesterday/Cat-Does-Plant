#deprecated code for making working files using custom scripts.

import pandas as pd
import os

######################## HQ ##########################


# dfa, dfl = pd.DataFrame(columns = ['CP']), pd.DataFrame(columns = ['CP'])

# for file in os.listdir('/home/ric/Cat-chan/SNPTmp2'):
    
#     df = pd.read_csv('/home/ric/Cat-chan/SNPTmp2/%s' %file)
    
#     keep = []
#     for i,r,a in df[['REF', 'ALT']].itertuples():
#         if ',' in a:
#             if len(r) == len(a.split(',')[0]) == len(a.split(',')[1]):
#                 keep.append(i)
#         else:
#             if len(r) == len(a):
#                 keep.append(i)
                
#     df = df.iloc[keep,:].reset_index(drop = True)
    
#     tdf1 = df[['CP', 'ALT']].copy()
#     tdf1 = tdf1.rename(columns = {'ALT' : file.split('.')[0]})
    
#     tdf2 = df[['CP', 'Alleles']].copy()
#     tdf2 = tdf2.rename(columns = {'Alleles' : file.split('.')[0]})
    
#     dfa = pd.merge(dfa, tdf1, on = 'CP', how = 'outer')
#     dfl = pd.merge(dfl, tdf2, on = 'CP', how = 'outer')
    
#     print('%d remain.' %(len(os.listdir('/home/ric/Cat-chan/SNPTmp2')) - os.listdir('/home/ric/Cat-chan/SNPTmp2').index(file)))
    
# dfa.to_csv('/home/ric/Cat-chan/WF1hqSNPa.csv', index = False) 
# dfl.to_csv('/home/ric/Cat-chan/WF1hqSNPl.csv', index = False)

######################## LQ ##########################


# dfa = pd.read_csv('/home/ric/Cat-chan/WF1hqSNPa.csv', dtype = str)[['CP']]
# dfl = pd.read_csv('/home/ric/Cat-chan/WF1hqSNPl.csv', dtype = str)[['CP']]

# for file in os.listdir('/home/ric/Cat-chan/BSNPrune'):
    
#     df = pd.read_csv('/home/ric/Cat-chan/BSNPrune/%s' %file)
    
#     tdf1 = df[['CP', 'ALT']].copy()
#     tdf1 = tdf1.rename(columns = {'ALT' : file.split('.')[0]})
    
#     tdf2 = df[['CP', 'Alleles']].copy()
#     tdf2 = tdf2.rename(columns = {'Alleles' : file.split('.')[0]})
    
#     dfa = pd.merge(dfa, tdf1, on = 'CP', how = 'left')
#     dfl = pd.merge(dfl, tdf2, on = 'CP', how = 'left')
    
#     print('%d remain.' %(len(os.listdir('/home/ric/Cat-chan/BSNPrune')) - os.listdir('/home/ric/Cat-chan/BSNPrune').index(file)))
    
# dfa.to_csv('/home/ric/Cat-chan/WF1lqSNPa.csv', index = False) 
# dfl.to_csv('/home/ric/Cat-chan/WF1lqSNPl.csv', index = False) 


####################### superimpose ####################

dfa = pd.read_csv('/home/ric/Cat-chan/WF1lqSNPa.csv', dtype = str)
for c in dfa.columns[1:]:
    dfa.loc[dfa[c].notnull(), c] = 'l'

pd.read_csv('/home/ric/Cat-chan/WF1hqSNPa.csv', dtype = str).combine_first(dfa).to_csv('/home/ric/Cat-chan/WF1hqlqSNPa.csv', index = False)

dfl = pd.read_csv('/home/ric/Cat-chan/WF1lqSNPl.csv', dtype = str)
for c in dfl.columns[1:]:
    dfl.loc[dfl[c].notnull(), c] = 'l'

pd.read_csv('/home/ric/Cat-chan/WF1hqSNPl.csv', dtype = str).combine_first(dfl).to_csv('/home/ric/Cat-chan/WF1hqlqSNPl.csv', index = False)

########################### tests #########################

df = pd.read_csv('WF1hqlqSNPa.csv', dtype = str)
df3 = pd.read_csv('WF1hqSNPa.csv', dtype = str)

assert sum(sum(df2.isnull().values)) == sum(sum((df == 'l').values) + sum((df.isnull()).values))