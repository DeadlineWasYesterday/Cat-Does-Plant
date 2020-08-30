import pandas as pd
from multiprocessing import Pool
import gzip, io, os, time

def read_vcf(path):
    with gzip.open(path, 'r') as f:
        lines = []
        for line in f:
            if line[:2] != b'##':
                l = line.decode('utf-8')
                if l.split('\t')[4] != '.' and l.split('\t')[6] != 'LowQual':
                    lines.append(l)

    lines = lines[1:]
    lines = ''.join(lines)
    return pd.read_csv(io.StringIO(lines), names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'FILE'],
           dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str,
                 'FILE': str}, sep='\t')


def poolparty(it):
    
    ft1 = time.time()
    
    df = read_vcf('/home/ric/Cat-chan/174vcf/%s' %it)

    ft2 = time.time() - ft1
    #print('Took %d seconds to read dataframe of shape %dx%d.' %((ft2,) + df.shape))       
    
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

    
    df = df.drop_duplicates(subset = 'CP').reset_index(drop = True)
    print('Shape before dropping %dx%d. Shape after dropping %dx%d.' % (bshape + df.shape))

    df.to_csv('/home/ric/Cat-chan/SNPTmp2/%s.csv' %it, index = False)
    #print('Write complete for %s.' %it)
    
    return 0


if __name__ == '__main__':

    t1 = time.time()
    print('Entered main function.')
    
    it = iter([f for f in os.listdir('/home/ric/Cat-chan/174vcf/') if f[-6:] == 'vcf.gz'])
        
    pool = Pool(processes = 15, maxtasksperchild = 20)    
    
    result = pool.map(poolparty, it)
    pool.close()    
    print('Pool closed.')    
   
    #rdf = pd.concat(result, axis = 0).reset_index(drop = True)      

    #save cdf
    #rdf.to_csv('/home/ric/Cat-chan/exonrdf.csv', index = False)
