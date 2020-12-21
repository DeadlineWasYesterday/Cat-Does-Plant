import pandas as pd
import numpy as np
import io


def func(fn):
    
    df = pd.read_csv('/home/ric/Cat-chan/SNPTmp2/%s' %fn)
    
    cl, pl, rl, al, ll = [],[],[],[],[]

    for q,r,a,l in df[['CP', 'REF', 'ALT', 'Alleles']].itertuples(index = False):

        if ',' in a:
            a1,a2 = a.split(',')
            if len(r) == len(a1) == len(a2) == 1:
                cl.append(q[:5])
                pl.append(q[5:])
                rl.append(r)
                al.append(a)
                ll.append(l)

            else: #indel
                la = max(a.split(','), key=len)
                if la == a1:
                    sa = a2
                else:
                    sa = a1            

                if len(r) > len(la): #del
                    la = la + 'N'
                    sa = sa + 'N'

                    for i in range(len(r)):
                        cl.append(q[:5])
                        pl.append(str(int(q[5:].split('r')[0]) + i))
                        rl.append(r[i])
                        al.append(sa[min(i, len(sa) - 1)] + ',' + la[min(i, len(la) - 1)])
                        ll.append(l) 

                else: #ins
                    nr = r + 'N'
                    sa = sa + 'N'
                    p = float(q[5:].split('i')[0])
                    for i in range(len(la)):
                        cl.append(q[:5])
                        pl.append(str(round(p,2)))
                        rl.append(nr[min(i, len(nr) - 1)])
                        al.append(sa[min(i, len(sa) - 1)] + ',' + la[min(i, len(la) - 1)])
                        ll.append(l)
                        p = p + 0.1

        else:        
            if len(r) == len(a) == 1:
                cl.append(q[:5])
                pl.append(q[5:])
                rl.append(r)
                al.append(a)
                ll.append(l)

            else: #indel              
                if len(r) >= len(a): #del
                    na = a + 'N'                    
                    for i in range(len(r)):
                        cl.append(q[:5])
                        pl.append(str(int(q[5:].split('r')[0]) + i))
                        rl.append(r[i])
                        al.append(na[min(i, len(na) - 1)])
                        ll.append(l)                                                                       

                else: #ins                    
                    nr = r + 'N'
                    p = float(q[5:].split('i')[0])
                    for i in range(len(a)):                    
                        cl.append(q[:5])
                        pl.append(str(round(p,2)))
                        rl.append(nr[min(i, len(nr) - 1)])
                        al.append(a[i])
                        ll.append(l)
                        p = p + 0.1


    df = pd.DataFrame(zip(cl,pl,rl,al,ll), columns = ['#CHROM', 'POS', 'REF', 'ALT', 'Alleles'])
    df['POS'] = df['POS'].astype(float)
    df['CP'] = df['#CHROM'] + df['POS'].astype(str)
    
    #fix r = a1/a2
    #fix a1 = a2
    
    kl,al,ll = [],[],[]
    for i,r,a,l in df[['REF','ALT','Alleles']].itertuples():
        if ',' in a:
            a1,a2 = a.split(',')
            
            if r == a1 == a2:
                kl.append(i)
                al.append(a)
                ll.append(l)                
            elif a1 == a2:
                al.append(a1)
                ll.append('1/1')
            elif r == a1:
                al.append(a2)
                ll.append('0/1')
            elif r == a2:
                al.append(a1)
                ll.append('0/1')                
            else:
                al.append(a)
                ll.append(l)
                
        else:
            if r == a:
                kl.append(i)
                al.append(a)
                ll.append(l)
                
            else:
                al.append(a)
                ll.append(l)
                
    df['ALT'], df['Alleles'] = al,ll

    df = df.iloc[df.index.difference(kl),:].reset_index(drop = True)

    df = df.sort_values('Alleles', ascending = False).drop_duplicates(['REF','ALT','CP']).reset_index(drop = True)

    tdf1 = df[df.duplicated('CP', keep = 'first')]
    tdf2 = df[df.duplicated('CP', keep = 'last')]

    df.loc[df.duplicated('CP'), 'ALT'] = [''.join(t) for t in zip(tdf1['ALT'], [','] * len(tdf1), tdf2['ALT'])]
    df.loc[df.duplicated('CP'), 'Alleles'] = '1/2'

    df = df.drop_duplicates('CP').reset_index(drop = True)