#Uses the merged VCF file from the 3000-rice genome project
# and makes a number of working files.
# WF1 = File with all SNPs with indels.
# WF2 = ALl biallelic SNPs with no indels.
# WF3 = Biallelic and Multiallelic Substitutions
# WF4 = Biallelic Substitutions and Indels
# WF5 = All substitutes and indels having phred-scaled quality > 30

#includes missing value (low-quality reads) filter and dealing with 
# an unique 'exclusive-1/1' SNP problem. 

#includes code for pruning files and keeping only 
#markers that are in approximate linkage equilibrium with one another.


import pandas as pd
import numpy as np
import io, os, bgzip, gzip, pickle


#From the merged VCF all the way upto file for input into Beagle 5.1.
#Not computationally efficient. Suitable for computers with at least 256GB of RAM.

with gzip.open('176merged.vcf.gz', 'r') as f:
    lines = []
    for line in f:
        if line[:2] != b'##':
            l = line.decode('utf-8')
            if l.split('\t')[4] != '.':
                lines.append(l)
#lines = list(map(lambda x: '\t'.join(list(map(lambda y: y.split(':')[0], x.split('\t')))), lines))


pickle.dump(lines, open("lines.pickle", "wb"))

lines = pickle.load(open('lines.pickle' , 'rb'))

df = pd.read_csv(io.StringIO(''.join(lines[:8000000])), dtype= str, sep='\t')

df['QUAL'] = '.'
df['INFO'] = '.'
df['FORMAT'] = '.'
df['ID'] = df.index
df = df.applymap(lambda x: x.split(':')[0])

df.to_csv('WF1p1.vcf', index = False, sep = '\t')


lines = pickle.load(open('lines.pickle' , 'rb'))

df = pd.read_csv(io.StringIO(''.join(lines[8000000:])), dtype= str, sep='\t', header = None, names = lines[0][:-1].split('\t'))

df['QUAL'] = '.'
df['INFO'] = '.'
df['FORMAT'] = '.'
df['ID'] = df.index + (8000000 - 1)
df['ID'] = df['ID'].astype(str)
df = df.applymap(lambda x: x.split(':')[0])

df.to_csv('WF1p2.vcf', index = False, sep = '\t')

df1 = pd.read_csv('WF1p1.vcf', dtype = str, sep = '\t')
df2 = pd.read_csv('WF1p2.vcf', dtype = str, sep = '\t')

df = pd.concat([df1,df2], axis = 0).reset_index(drop = True)
df.to_csv('WF1BIn.vcf', index = False, sep = '\t')

#Beagle input file WF1BIn.vcf has been written



#Beagle output file is is WF1.vcf.gz
#We are converting phased VCF to unphased since 3KRG data is unphased.
with gzip.open('WF1.vcf.gz', 'r') as f:
    lines = []
    for line in f:
        if line[:2] != b'##':
            lines.append(line.decode('utf-8').replace('|', '/'))
df = pd.read_csv(io.StringIO(''.join(lines)), dtype= str, sep='\t')

df.to_csv('WF1BOut.vcf', index = False, sep = '\t')
#Unphased VCF file has been written.


#Filtering alleles having missing rates greater than 20%.
#The 20% mark has been hardcoded here by 70 missing alleles.
df2 = pd.read_csv('WF1BIn.vcf', sep = '\t')

keep = []
for i,r in df2.iterrows():
    s = (''.join(r[-176:].to_list()))
    if (s.count('.')) <= 70:
        keep.append(i)
        
df = df.loc[keep, :]
df.to_csv('WF1M20.vcf', sep = '\t', index = False)
#File with SNPs with >20% missing values pruned out has been written. 


### Reference allele filter

#This step is necessary for the removal of SNPs that have no or very few 0/0 alleles.
#SNPs with minimal number of 0/0 genotypes to keep has been hardcoded as 8. Much like a 5% MAF filter.
df = pd.read_csv('WF1M20.vcf', sep = '\t', dtype = str)
keep = []
for i,r in df.iterrows():
    if ''.join(r).count('0/0') > 8:
        keep.append(i)
        
df = df.loc[keep, :]
df.to_csv('WF1M20R8.vcf', index = False, sep = '\t')
#File with 'exclusive-1/1' issue fixed has been written.



### WF2 Biallelic Substitutions Only
df = pd.read_csv('WF1M20R8.vcf', sep = '\t', dtype = str)

df['t'] = df['REF'].apply(len) + df['ALT'].apply(len)
df = df[df['t'] == 2].reset_index(drop = True)
df = df.drop(columns = ['t'])

df.to_csv('WF2.vcf', sep = '\t', index = False)

### WF3 Biallelic and Multiallelic Substitutions
df = pd.read_csv('WF1M20R8.vcf', sep = '\t', dtype = str)

df['t'] = df['ALT'].apply(lambda x: len(max(x.split(','), key = len)))
df['t'] = df['REF'].apply(len) + df['t']
df = df[df['t'] == 2].reset_index(drop = True)
df = df.drop(columns = ['t'])

df.to_csv('WF3.vcf', sep = '\t', index = False)

### WF4 Biallelic Substitutions and Indels
df = pd.read_csv('WF1M20R8.vcf', sep = '\t', dtype = str)

df['t'] = df['ALT'].apply(lambda x: len(x.split(',')))
df = df[df['t'] == 1]
df = df.drop(columns = ['t'])

df.to_csv('WF4.vcf', sep = '\t', index = False)


### WF5 HQ Subs and Indels only
df = pd.read_csv('WF1BOut.vcf', dtype = str, sep = '\t')
df2 = pd.read_csv('WF1BIn.vcf', dtype = str, sep = '\t')


keep = []
for i,r in df2.iterrows():
    s = (''.join(r[-176:].to_list()))
    if '.' not in s:
        keep.append(i)
        
df = df.loc[keep, :]

df.to_csv('WF5.vcf', sep = '\t', index = False)
#ref allele filter
keep = []
for i,r in df.iterrows():
    if ''.join(r).count('0/0') > 8:
        keep.append(i)
        
df = df.loc[keep, :]
df.to_csv('WF5R8.vcf', index = False, sep = '\t')




###Pruned WF2

#The prune.in files are from PLINK. Refer to shell scripts.

df = pd.read_csv('WF2.vcf', dtype = str, sep = '\t')
df.index = df['ID'].values.astype(int)
df = df.loc[pd.read_csv('WF2.prune.in', header = None)[0].to_list(), :]

df.to_csv('WF2LE.vcf', index = False, sep = '\t')


###Pruned WF5R8

df = pd.read_csv('Data/WF5R8.vcf', dtype = str, sep = '\t')
df.index = df['ID'].astype(int).values
df.loc[pd.read_csv('Data/WF5R8.prune.in', header = None)[0],:].to_csv('Data/WF5R8LE.vcf', sep = '\t', index = False)



###Output numeric
#Numeric output facilitates the process of making genotype matrices

df = pd.read_csv('Data/WF2.vcf', dtype = str, sep = '\t', header = None)
df2 = df.applymap(lambda x: x.count('1'))
df2 = df2.iloc[:, 8:]
df2.iloc[:,0] = df.loc[:, 2]

df2.iloc[0,:] = ['Taxa', 'IRIS_313-10340',
 'IRIS_313-10341',
 'IRIS_313-10537',
 'IRIS_313-10539',
 'IRIS_313-10587',
 'IRIS_313-10592',
 'IRIS_313-10593',
 'IRIS_313-10595',
 'IRIS_313-10598',
 'IRIS_313-10600',
 'IRIS_313-10602',
 'IRIS_313-10603',
 'IRIS_313-10604',
 'IRIS_313-10605',
 'IRIS_313-10606',
 'IRIS_313-10963',
 'IRIS_313-10964',
 'IRIS_313-10965',
 'IRIS_313-10971',
 'IRIS_313-10972',
 'IRIS_313-10973',
 'IRIS_313-10974',
 'IRIS_313-10975',
 'IRIS_313-10976',
 'IRIS_313-10977',
 'IRIS_313-10978',
 'IRIS_313-10979',
 'IRIS_313-10980',
 'IRIS_313-10981',
 'IRIS_313-10982',
 'IRIS_313-10983',
 'IRIS_313-10985',
 'IRIS_313-10986',
 'IRIS_313-10987',
 'IRIS_313-11013',
 'IRIS_313-11014',
 'IRIS_313-11015',
 'IRIS_313-11016',
 'IRIS_313-11017',
 'IRIS_313-11018',
 'IRIS_313-11019',
 'IRIS_313-11047',
 'IRIS_313-11048',
 'IRIS_313-11049',
 'IRIS_313-11050',
 'IRIS_313-11052',
 'IRIS_313-11053',
 'IRIS_313-11054',
 'IRIS_313-11056',
 'IRIS_313-11057',
 'IRIS_313-11058',
 'IRIS_313-11059',
 'IRIS_313-11060',
 'IRIS_313-11061',
 'IRIS_313-11062',
 'IRIS_313-11063',
 'IRIS_313-11064',
 'IRIS_313-11065',
 'IRIS_313-11066',
 'IRIS_313-11067',
 'IRIS_313-11068',
 'IRIS_313-11069',
 'IRIS_313-11070',
 'IRIS_313-11111',
 'IRIS_313-11112',
 'IRIS_313-11113',
 'IRIS_313-11114',
 'IRIS_313-11115',
 'IRIS_313-11116',
 'IRIS_313-11123',
 'IRIS_313-11124',
 'IRIS_313-11154',
 'IRIS_313-11163',
 'IRIS_313-11203',
 'IRIS_313-11204',
 'IRIS_313-11205',
 'IRIS_313-11206',
 'IRIS_313-11207',
 'IRIS_313-11208',
 'IRIS_313-11210',
 'IRIS_313-11211',
 'IRIS_313-11212',
 'IRIS_313-11213',
 'IRIS_313-11214',
 'IRIS_313-11216',
 'IRIS_313-11217',
 'IRIS_313-11218',
 'IRIS_313-11219',
 'IRIS_313-11220',
 'IRIS_313-11221',
 'IRIS_313-11222',
 'IRIS_313-11223',
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
 'IRIS_313-11322',
 'IRIS_313-11323',
 'IRIS_313-11324',
 'IRIS_313-11326',
 'IRIS_313-11399',
 'IRIS_313-11400',
 'IRIS_313-11401',
 'IRIS_313-11402',
 'IRIS_313-11404',
 'IRIS_313-11481',
 'IRIS_313-11482',
 'IRIS_313-11483',
 'IRIS_313-11484',
 'IRIS_313-11486',
 'IRIS_313-11487',
 'IRIS_313-11557',
 'IRIS_313-11558',
 'IRIS_313-11722',
 'IRIS_313-11945',
 'IRIS_313-11964',
 'IRIS_313-12002',
 'IRIS_313-12055',
 'IRIS_313-12094',
 'IRIS_313-12141',
 'IRIS_313-8244',
 'IRIS_313-8252',
 'IRIS_313-8283',
 'IRIS_313-8321',
 'IRIS_313-8349',
 'IRIS_313-8410',
 'IRIS_313-8437',
 'IRIS_313-8509',
 'IRIS_313-8632',
 'IRIS_313-8641',
 'IRIS_313-8683',
 'IRIS_313-8703',
 'IRIS_313-8717',
 'IRIS_313-8721',
 'IRIS_313-8733',
 'IRIS_313-8737',
 'IRIS_313-8789',
 'IRIS_313-8822',
 'IRIS_313-8850',
 'IRIS_313-8854',
 'IRIS_313-8864',
 'IRIS_313-8930',
 'IRIS_313-8932',
 'IRIS_313-8963',
 'IRIS_313-9049',
 'IRIS_313-9066',
 'IRIS_313-9067',
 'IRIS_313-9072',
 'IRIS_313-9108',
 'IRIS_313-9139',
 'IRIS_313-9148',
 'IRIS_313-9156',
 'IRIS_313-9170',
 'IRIS_313-9174',
 'IRIS_313-9218',
 'IRIS_313-9262',
 'IRIS_313-9271',
 'IRIS_313-9325',
 'IRIS_313-9368',
 'IRIS_313-9391',
 'IRIS_313-9397',
 'IRIS_313-9422',
 'IRIS_313-9445',
 'IRIS_313-9594',
 'IRIS_313-9606',
 'IRIS_313-9617',
 'IRIS_313-9626',
 'IRIS_313-9661',
 'IRIS_313-9682',
 'IRIS_313-9696']

np.savetxt('WF2.num.txt', df2.transpose().values, fmt = '%s', delimiter="\t")