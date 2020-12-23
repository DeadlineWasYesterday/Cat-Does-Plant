#Tries to compile a vcf file from a dataframe of substitutions and a dataframe of genotype (./.)

#Has bugs.
#Three types of exceptions
#Cells with commas
#SNPs leading into insertions
#SNPs leading into deletions

import os, pickle
import pandas as pd
from multiprocessing import Pool


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
        

def efunc(fname):
    
    f = open('/home/ric/Cat-chan/rcvtin/%s' %fname, 'rb')
    it = pickle.load(f)
    f.close()
    
    print('%s loaded. Shape %dx%d, %dx%d.' % ((fname,) + it[0].shape + it[1].shape))
    
    rdf = pd.DataFrame(columns = list(it[0].columns) + ['ALT'])
    
    it = (it[0].reset_index(drop = True), it[1].reset_index(drop = True))

    for i in it[0].index:

        s1 = it[0].iloc[i,:][1:].copy()
        s2 = it[1].iloc[i,:][1:].copy()

        als = []
        for sub in s1.dropna().unique():
            if ',' in sub:
                als = als + sub.split(',')
            else:
                als.append(sub)
        als = list(set(als))

        for sub in s1.dropna().unique():
            for p, g in s2[s1 == sub].iteritems():
                if g.split('/')[0] != '0':
                    t1 = als.index(sub.split(',')[int(g.split('/')[0]) - 1])
                else:
                    t1 = 0

                if g.split('/')[1] != '0':
                    t2 = als.index(sub.split(',')[int(g.split('/')[1]) - 1])
                else:
                    t2 = 0

                s2.loc[p] = str(t1 + 1) + '/' + str(t2 + 1)

        s2 = it[0].iloc[i,:][:1].append(s2).append(pd.Series({'ALT' : ','.join(als)}))        
        s2.name = i
        rdf = rdf.append(s2)
              
        
    rdf.to_csv('/home/ric/Cat-chan/rvctemp/%s.csv' %fname, index = False)  
    print('Write complete for %s.' % fname)
    
    return rdf

if __name__ == '__main__':
    
    df1 = pd.read_csv('/home/ric/Cat-chan/WF1hqSNPa.csv', dtype = str)
    df2 = pd.read_csv('/home/ric/Cat-chan/WF2hqSNPl.csv', dtype = str)
    
    lists = list(chunks(df1.index.to_list(), 30000))

    for lst in lists:
        f = open('/home/ric/Cat-chan/rcvtin/%s.pkl' %lists.index(lst), 'wb')
        pickle.dump((df1.iloc[lst,:].copy(),df2.iloc[lst,:].copy()), f)
        f.close()

    print('Initiating pool.')
    print('%d files.' % len(os.listdir('/home/ric/Cat-chan/rcvtin')))
    pool = Pool(processes = 15, maxtasksperchild = 30)
    
    rs = pool.map(efunc, os.listdir('/home/ric/Cat-chan/rcvtin'))
    
    pool.close()
    
    print('Pool closed.')
    
    pd.concat(rs, axis = 0).to_csv('/home/ric/Cat-chan/WF1RVC.vcf', index = False)
