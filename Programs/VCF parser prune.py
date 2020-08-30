import pandas as pd
from multiprocessing import Pool
import gzip, io, os, time

def read_vcf(path):
    with gzip.open(path, 'r') as f:
        lines = []
        for line in f:
            if line[:2] != b'##':
                l = line.decode('utf-8')
                if l.split('\t')[6] == 'LowQual':
                    lines.append(l)

    lines = ''.join(lines)
    return pd.read_csv(io.StringIO(lines), names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'FILE'],
           dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str,
                 'FILE': str}, sep='\t')

def poolparty(tit):
    
    it, pdf = tit    
    
    ft1 = time.time()

    df = read_vcf('/home/ric/Cat-chan/174vcf/%s' %it)

    ft2 = time.time() - ft1
    print('Took %d seconds to read dataframe of shape %dx%d.' %((ft2,) + df.shape))       
    
    df['Alleles'] = [s.split(':')[0] for s in df.iloc[:, -1].to_list()]
    df['CP'] = df['#CHROM'] + df['POS'].astype(str)
    df = df[['CP', 'REF', 'ALT', 'Alleles']]  

    bshape = df.shape
    
    for pos in df['CP'][df['CP'].duplicated()].to_list():
        if len(df.loc[df['CP'] == pos, 'REF'].iloc[0]) + len(df.loc[df['CP'] == pos, 'REF'].iloc[1]) == 2:
            df.loc[df['CP'] == pos, 'CP'] = [pos, pos + 'i']
        else:
            #print('Del ALTs: %s.' % df.loc[df['CP'] == pos, 'ALT'].iloc[0] + df.loc[df['CP'] == pos, 'ALT'].iloc[1])
            df.loc[df['CP'] == pos, 'CP'] = [pos, pos + 'd']

    
    df = df.drop(columns = ['REF'])
    df = df.drop_duplicates(subset = 'CP').reset_index(drop = True)
    print('Shape before dropping %dx%d. Shape after dropping %dx%d.' % (bshape + df.shape))
    
    pdf = pd.merge(pdf, df, on ='CP', how = 'left').reset_index(drop = True)

    tf3 = time.time() - (ft1 + ft2)
    print('Merge complete for %s in %d seconds. Shape: %dx%d.' %(it, tf3, pdf.shape[0], pdf.shape[1]))
            
    pdf.to_csv('/home/ric/Cat-chan/BSNPrune/%s' %it, index = False) 
      

if __name__ == '__main__':

    t1 = time.time()
    print('Entered working loop.')
    
    pdf = pd.read_csv('/home/ric/Cat-chan/hqBSNP.csv', dtype = str)[['CP']]
    
    tit = iter([(f,pdf) for f in os.listdir('/home/ric/Cat-chan/174vcf/') if f[-6:] == 'vcf.gz'])
    
    pool = Pool(processes = 15, maxtasksperchild = 20)
    
    result = pool.map(poolparty, tit)

    pool.close()
    
    print('Pool closed.')
