import pandas as pd
import numpy as np
from multiprocessing import Pool

def func(it):
    n = np.array([], dtype = 'int')
    df = pd.read_csv('pan vs exon 100, chronly.csv')
    df = df[df['1'] == it]
    
    print('Inside funcit. Data loaded. Working with %dx%d.' %df.shape) 
    
    for s,e in df[['8','9']].itertuples(index = False):
        if s < e:
            t = np.arange(s, e+1, dtype = 'int')
            n = np.concatenate([n, t[~np.isin(t,n)]])
        else:
            t = np.arange(e, s+1, dtype = 'int')
            n = np.concatenate([n, t[~np.isin(t,n)]])

    print('Done with %s.' %it)
    
    return list(zip([it] * len(n), list(n)))
    


if __name__ == '__main__':

    n = np.array([], dtype = 'int')
    
    it = iter(['chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07',
       'chr08', 'chr09', 'chr10', 'chr11', 'chr12'])
    
    pool = Pool(processes = 12)
    
    result = pool.map(func, it)
    
    print('Result generated.')
    
    tlst = np.concatenate(result)
                
    rdf = pd.DataFrame(tlst, columns = ['Chromosome', 'Position'])
                
    rdf.to_csv('RAPDB exon position map.csv', index = False)
                
