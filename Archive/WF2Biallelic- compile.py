#deprecated code for making working files using custom scripts.

import pandas as pd
import numpy as np
import os

################### compile alleles ########################

# dfa, dfl = pd.DataFrame(columns = ['CP']), pd.DataFrame(columns = ['CPRA'])

# for file in os.listdir('/home/ric/Cat-chan/SNPTmp2'):
    
#     df = pd.read_csv('/home/ric/Cat-chan/SNPTmp2/%s' %file)
    
#     df['CPRA'] = df['CP'] + '-' + df['REF'] + '-' + df['ALT']
    
#     tdf1 = df[['CPRA', 'Alleles']].copy()
#     tdf1 = tdf1.rename(columns = {'Alleles' : file.split('.')[0]})
    
#     dfl = pd.merge(dfl, tdf1, on = 'CPRA', how = 'outer')
    
#     print('%d remain.' %(len(os.listdir('/home/ric/Cat-chan/SNPTmp2')) - os.listdir('/home/ric/Cat-chan/SNPTmp2').index(file)))
    
# dfl.to_csv('/home/ric/Cat-chan/WF2Biallelicl.csv', index = False) 


##################### drop multiallelic ###################


# dfl = pd.read_csv('/home/ric/Cat-chan/WF2Biallelicl.csv', dtype = str)

# dfl['CP'] = [s.split('-')[0] for s in dfl['CPRA'].to_list()]

# cp = []
# for c in dfl['CP'].to_list():
#     if c[5:].isdigit():
#         cp.append(c)
#     elif c[5:][:-1].isdigit():
#         cp.append(c[:-1])
#     else:
#         cp.append(c.split('r')[0])
        
# dfl['CP'] = cp

# dfl = dfl.iloc[dfl.index.difference(dfl[dfl.duplicated('CP', keep = False)].index), :].reset_index(drop = True)

# dfl.to_csv('/home/ric/Cat-chan/WF2BialleliclID.csv', index = False) 


############################## Biallelic SNP only #################################

# dfl = pd.read_csv('/home/ric/Cat-chan/WF2BiallelicID.csv', dtype = str)

# dfl['REF'] = [s.split('-')[1] for s in dfl['CPRA'].to_list()]
# dfl['ALT'] = [s.split('-')[2] for s in dfl['CPRA'].to_list()]

# keep = []
# for i,r,a in dfl[['REF', 'ALT']].itertuples():
#     if len(r) == len(a) == 1:
#         keep.append(i)
        
# dfl = dfl.iloc[keep, :].reset_index(drop = True)

# #drop ubiquitous
# keep = []
# for i,r in dfl.iterrows():
#     if len(r[1:-3].unique()) > 1:
#         keep.append(i)
    
# dfl = dfl.iloc[keep, :]

# dfl.to_csv('/home/ric/Cat-chan/WF2BialleliclSNP.csv', index = False) 



######################### method 2: from WF1 #########################


# dfa = pd.read_csv('/home/ric/Cat-chan/WF1hqSNPa.csv', dtype = str)
# dfl = pd.read_csv('/home/ric/Cat-chan/WF1hqSNPl.csv', dtype = str)

# keep = []
# for i,r in dfa.iterrows():
#     if len(r.unique()) == 3:
#         keep.append(i)

# dfa = dfa.loc[keep,:].reset_index(drop = True)
# dfl = dfl.loc[keep,:].reset_index(drop = True)

# keep = []
# for i,r in dfl.iterrows():
#     if '1/2' not in r.unique():
#         keep.append(i)
        
# dfa = dfa.loc[keep,:].reset_index(drop = True)
# dfl = dfl.loc[keep,:].reset_index(drop = True)

# alt = []
# for i,r in dfa.iterrows():
#     alt.append(r.iloc[1:].dropna().unique()[0])
    
# dfl['ALT'] = alt

# dfl.to_csv('/home/ric/Cat-chan/WF2BiallelicM2.csv', index = False)


############################ endow LQ cells ######################

dfl = pd.read_csv('/home/ric/Cat-chan/WF2BialleliclSNP.csv', dtype = str)
df = dfl[['CP']].copy()

for file in os.listdir('/home/ric/Cat-chan/BSNPrune/'):
    tdf = pd.read_csv('/home/ric/Cat-chan/BSNPrune/%s' %file)[['CP', 'ALT']]
    tdf.loc[tdf['ALT'].notnull(), 'ALT'] = 'l'
    tdf = tdf.rename(columns = {'ALT' : file.split('.')[0]})
    
    df = pd.merge(df, tdf, on = 'CP', how = 'left')
    
dfl = dfl.combine_first(df)

dfl.to_csv('/home/ric/Cat-chan/WF2BialleliclSNPwlq.csv', index = False)



############################ superimpose test ###############################

df = pd.read_csv('WF2BialleliclSNP.csv', dtype = str)
df2 = pd.read_csv('WF2BialleliclSNPwlq.csv', dtype = str)

assert sum(sum(df.isnull().values)) == sum(sum((df2 == 'l').values) + sum((df2.isnull()).values))



############################# to VCF ######################################

df = pd.read_csv('WF2BialleliclSNPwlq.csv', dtype = str)

df = df.fillna('0/0')

df[df == 'l'] = np.nan

df['#CHROM'] = [s[:5] for s in df['CP'].to_list()]
df['POS'] = [int(s[5:]) for s in df['CP'].to_list()]

df['ID'] = df.index
df['QUAL'] = '.'
df['FILTER'] = '.'
df['INFO'] = '.'
df['FORMAT'] = 'GT'

df.loc[:, ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + list(df.columns[3:-8])].to_csv('WF2Bial.vcf', index = False, sep = '\t')


############################## imputation test #################################

df = pd.read_csv('WF2Bial.vcf', sep = '\t')
dfo = df.copy()

print(sum(sum(df.isnull().values)))
#108744649 null values

print(108744649 / 174)
#624969.2471264368 column mean

df = df.dropna(axis = 0).reset_index(drop = True)
#3.4M hq rows. Originally 13.35M

print(625000 * 3.5 / 13.4)
#163246.26865671642 null values expected per hq only column


for c in df.columns[9:]:
    df.loc[set(np.random.randint(0, 3400831, 175000)), c] = np.nan
    
print(sum(sum(df.isnull().values)))
#29679143 total null values generated 

df = df.fillna('./.')

df.to_csv('WF2ActestIn.vcf', index = False, sep = '\t')

#outfile read function
def read_vcf(path):
    with open(path, 'r') as f:
        lines = []
        for line in f:
            if line[:2] != '##':
                lines.append(line.replace('|', '/'))


    lines = ''.join(lines)
    return pd.read_csv(io.StringIO(lines), sep='\t')


df = pd.read_csv('WF2Bial.vcf', sep = '\t')
df = df.dropna(axis = 0).reset_index(drop = True)


######with defaults
df2 = read_vcf('WF2outdef.vcf')
df2 = df2.replace('1/0', '0/1')
df2['FILTER'] = '.'

print(sum(sum((df != df2).values)))
#255383

######with ne=500000
df2 = read_vcf('WF2outne5.vcf')
df2 = df2.replace('1/0', '0/1')
df2['FILTER'] = '.'

print(sum(sum((df != df2).values)))
#257566



######with ne=750000
df2 = read_vcf('WF2outne75.vcf')
df2 = df2.replace('1/0', '0/1')
df2['FILTER'] = '.'

print(sum(sum((df != df2).values)))
#256231

######with ne=2000000
df2 = read_vcf('WF2outne2.vcf')
df2 = df2.replace('1/0', '0/1')
df2['FILTER'] = '.'

print(sum(sum((df != df2).values)))
#254176







##################### Beagle input file #################

df = pd.read_csv('WF2Bial.vcf', sep = '\t').fillna('./.')

df['POS'] = df['POS'].astype(int)
df = df.sort_values(by = ['#CHROM', 'POS']).reset_index(drop = True)
df['ID'] = df.index

df.to_csv('WF2BialBeagIn.vcf', index = False, sep = '\t')






########################### Order Beagle Out ###########################






############### hmp file from beagle out ############

df = read_vcf('WF2BeagOut.vcf')

t = []
for i,r in df.iterrows():
    r[r == '0/0'] = r['REF'] + r['REF']
    r[r == '0/1'] = r['REF'] + r['ALT']
    r[r == '1/0'] = r['ALT'] + r['REF']
    r[r == '1/1'] = r['ALT'] + r['ALT']
    t.append(r)
    
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

df2.to_csv('WF2.hmp.txt', sep = '\t', index = False)

df2 = pd.read_csv('WF2.hmp.txt', sep = '\t', header = None)
np.savetxt('WF2np.hmp.txt', df2.fillna('NA').values, fmt = '%s', delimiter="\t")  




##################### MAF filter ####################

#Min MAF 12
df = pd.read_csv('WF2.hmp.txt', sep = '\t', dtype = str)

keep = []
for i,r in df.iterrows():
    s = (''.join(r[-174:].to_list()))
    if 41 < (s.count(s[0])) < 307:
        keep.append(i)
        
df2 = df.iloc[keep, :]
df2.to_csv('WF24MhmpMAF12.csv', index = False)
df2 = pd.read_csv('WF24MhmpMAF12.csv', header = None, dtype = str)
np.savetxt('Data/WF2MAF12np.hmp.txt', df2.fillna('NA').values, fmt = '%s', delimiter="\t")  

#Min MAF 15

df = pd.read_csv('WF2.hmp.txt', sep = '\t', dtype = str)

keep = []
for i,r in df.iterrows():
    s = (''.join(r[-174:].to_list()))
    if 53 < (s.count(s[0])) < 295:
        keep.append(i)
        
df2 = df.iloc[keep, :]
df2.to_csv('WF24MhmpMAF15.csv', index = False)
df2 = pd.read_csv('WF2hmpMAF15.csv', header = None, dtype = str)
np.savetxt('Data/WF2MAF15np.hmp.txt', df2.fillna('NA').values, fmt = '%s', delimiter="\t")  

