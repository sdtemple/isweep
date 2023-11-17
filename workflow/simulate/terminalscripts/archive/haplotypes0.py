# haplotype frequency, position
# seth d temple, sdtemple@uw.edu
# june 17, 2023

import sys
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt

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
 scocol,
 numsnp)=sys.argv[1:]

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
               afsize=0.100,
               afstep=0.025,
               afcol='AAF',
               winsize=0.5,
               winstep=0.1,
               wincol='CM',
               scorecol='SCORE',
               numsnp=5
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
    def mean(x):
        return sum(x)/len(x)
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
                haplotable.sort_values(by=scorecol,ascending=False,inplace=True)
                score=mean(list(haplotable[scorecol])[:numsnp])
                freq=mean(list(haplotable[afcol])[:numsnp])
                posi=mean(list(haplotable[wincol])[:numsnp])
                tup=(winmid,
                     aafmid,
                     score,
                     haplotable.shape[0],
                     freq,
                     posi
                    )
                windowed[ctr]=list(tup)
                ctr+=1
    windowed=pd.DataFrame(windowed)
    windowed=windowed.T
    windowed.columns=['POS','AF','SCORE','SIZE','FREQ','LOC']
    windowed.sort_values(by='SCORE',ascending=False,inplace=True)
    windowed.reset_index(inplace=True,drop=True)
    return windowed

# haplotype math
snptab=pd.read_csv(snpin,sep='\t')
ibdtab=pd.read_csv(mapin,sep='\t')
snptabcm=centiMorgan(snptab,ibdtab[bpcol],ibdtab[cmcol],bpidx)
haptab=haplotypes(snptabcm,
                  afsize,
                  afstep,
                  afcol,
                  winsize,
                  winstep,
                  wincol,
                  scorecol=scocol
                 )
haptab.sort_values(by='LOC',ascending=True,inplace=True)
haptabbp=basepair(haptab,ibdtab[bpcol],ibdtab[cmcol],-1)
haptabbp.sort_values(by='SCORE',ascending=False,inplace=True)
haptabbp.to_csv(folderout+'/haplotypes.tsv',sep='\t',index=False)

# best haplotype
besthap=list(haptabbp.iloc[0])
bestbp=besthap[-1]
bestcm=besthap[-2]
bestaf=besthap[-3]

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

# haplotypes figure
plt.scatter(haptab['POS'],haptab['FREQ'],c=haptab['SCORE'],cmap='copper_r',s=haptab['SIZE'])
plt.ylim(-0.1,1.1)
plt.colorbar(label='Haplotype z-score')
plt.ylabel('Haplotype frequency')
plt.xlabel('Position')
plt.savefig(folderout+'/isweep.hap.png',dpi=300)
plt.clf()

# snps figure
plt.scatter(snptabcm[wincol],snptabcm[afcol],c=snptabcm[scocol],cmap='copper_r',s=5)
plt.ylim(-0.1,1.1)
plt.colorbar(label='SNP z-score')
plt.ylabel('SNP frequency')
plt.xlabel('Position')
plt.savefig(folderout+'/isweep.snp.png',dpi=300)
plt.clf()
