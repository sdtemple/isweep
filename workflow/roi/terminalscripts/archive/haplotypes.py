# haplotype frequency, position
# seth d temple, sdtemple@uw.edu
# june 17, 2023

import sys
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

(snpin,
 mapin,
 folderout,
 bpidx,
 bpcol,
 cmcol,
 afsize,
 afstep,
 afcol,
 winsize,
 winstep,
 wincol,
 scocol)=sys.argv[1:]

winsize=float(winsize)
winstep=float(winstep)
afsize=float(afsize)
afstep=float(afstep)
bpidx=int(float(bpidx))

def centiMorgan(table,
                bpcol,
                cmcol,
                idx=0
               ):
    '''Append a column for cM distance

    Parameters
    ----------
    table : pandas DataFrame
        Rows are SNPs with a bp column
    bpcol : array-like
        bp column
    cmcol : array-like
        cM column (match bpcol)
    idx : int
        Column index with bp in table

    Returns
    -------
    pandas DataFrame
        Same table but with a cM column
    '''
    nm=list(table.columns)[idx]
    z1=np.sort(bpcol)
    z2=np.sort(cmcol)
    X=[]
    Y=table[nm]
    for y in Y:
        bp1=z1[z1<=y].tolist()[-1]
        bp2=z1[z1>y].tolist()[0]
        cm1=z2[z1<=y].tolist()[-1]
        cm2=z2[z1>y].tolist()[0]
        c=bp1+bp2
        a=bp1/c
        b=bp2/c
        cm=cm1*a+cm2*b
        X.append(cm)
    table['CM']=X
    return table

def basepair(table,
             bpcol,
             cmcol,
             idx=0
            ):
    '''Append a column for bp distance

    Parameters
    ----------
    table : pandas DataFrame
        Rows are SNPs with a bp column
    bpcol : array-like
        bp column
    cmcol : array-like
        cM column (match bpcol)
    idx : int
        Column index with bp in table

    Returns
    -------
    pandas DataFrame
        Same table but with a cM column
    '''
    nm=list(table.columns)[idx]
    z1=np.sort(bpcol)
    z2=np.sort(cmcol)
    X=[]
    Y=np.sort(table[nm])
    for y in Y:
        cm1=z2[z2<=y].tolist()[-1]
        cm2=z2[z2>y].tolist()[0]
        bp1=z1[z2<=y].tolist()[-1]
        bp2=z1[z2>y].tolist()[0]
        c=cm1+cm2
        a=cm1/c
        b=cm2/c
        bp=bp1*a+bp2*b
        X.append(bp)
    table['BP']=X
    return table

def haplotypes(table,
               afsize=0.05,
               afstep=0.025,
               afcol='AAF',
               winsize=0.2,
               winstep=0.1,
               wincol='CM',
               scorecol='SCORE'
              ):
    '''Define haplotypes

    Parameters
    ----------

    Returns
    -------
    pandas.DataFrame
        Haplotype table
    '''
    # local function
    def double_window(table,
                      leftaaf,
                      rightaaf,
                      aafcol,
                      leftwin,
                      rightwin,
                      wincol
                     ):
        '''Grid subset of table (frequency by position)'''
        subtable=table[(table[aafcol]>=leftaaf)&(table[aafcol]<=rightaaf)]
        subtable=subtable[(table[wincol]>=leftwin)&(subtable[wincol]<=rightwin)]
        return subtable
    # code
    winmin=table[wincol].min()
    winmax=table[wincol].max()
    aafmin=0
    aafmax=1
    numsnp=10
    aafrange=np.arange(aafmin,aafmax,afstep)
    winrange=np.arange(winmin,winmax,winstep)
    windowed=dict()
    ctr=0
    for a in aafrange:
        aafleft=a
        aafright=a+afsize
        aafmid=(aafleft+aafright)/2
        for w in winrange:
            winleft=w
            winright=w+winsize
            winmid=(winleft+winright)/2
            haplotable=double_window(table,
                                     aafleft,
                                     aafright,
                                     afcol,
                                     winleft,
                                     winright,
                                     wincol
                                    )
            if haplotable.shape[0] >= numsnp:
                tup=(winmid,
                     aafmid,
                     haplotable[scorecol].mean(),
                     haplotable.shape[0]
                    )
                windowed[ctr]=list(tup)
                ctr+=1
    windowed=pd.DataFrame(windowed)
    windowed=windowed.T
    windowed.columns=['POS','AF','SCORE','SIZE']
    windowed.sort_values(by='SCORE',ascending=False,inplace=True)
    windowed.reset_index(inplace=True,drop=True)
    return windowed

# haplotype math
snptab=pd.read_csv(snpin,sep='\t')
ibdtab=pd.read_csv(mapin,sep='\t')
snptabcm=centiMorgan(snptab,ibdtab[bpcol],ibdtab[cmcol],bpidx)
haptab=haplotypes(snptabcm,scorecol=scocol)
haptab.sort_values(by='POS',ascending=True,inplace=True)
haptabcm=basepair(haptab,ibdtab['BPWINDOW'],ibdtab['CMWINDOW'],0)
haptabcm.sort_values(by='SCORE',ascending=False,inplace=True)
haptabcm.to_csv(folderout+'/isweep.haplotypes.tsv',sep='\t',index=False)

# best haplotype
besthap=list(haptabcm.iloc[0])
bestbp=besthap[-1]
bestcm=besthap[0]
bestaf=besthap[1]

# frequency
f=open(folderout+'/isweep.hap.freq.txt','w')
f.write(str(round(bestaf,4)))
f.write('\n')
f.close()

# position
f=open(folderout+'/isweep.hap.pos.txt','w')
f.write(str(int(bestbp)))
f.write('\t')
f.write(str(round(bestcm,4)))
f.write('\n')
f.close()
