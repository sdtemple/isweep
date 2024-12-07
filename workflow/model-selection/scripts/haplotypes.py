# estimate haplotype frequency, position

import sys
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt

(snpin,
 folderout,
 winidx,
 freqidx,
 scoreidx,
 freqsize,
 freqstep,
 winsize,
 winstep,
 numsnp,
 low
 )=sys.argv[1:]

winsize=float(winsize)
winstep=float(winstep)
freqsize=float(freqsize)
freqstep=float(freqstep)

winidx=int(float(winidx))
freqidx=int(float(freqidx))
scoreidx=int(float(scoreidx))

low=float(low)
numsnp=int(float(numsnp))

def haplotypes(table,
               freqsize=0.04,
               freqstep=0.02,
               freqidx=1,
               winsize=200_000,
               winstep=50_000,
               winidx=0,
               scoreidx=-1,
               numsnp=5
              ):
    '''Define haplotypes

    Parameters
    ----------
    table : pandas.DataFrame
    freqsize : float
        Size of frequency window
    freqstep : float
        Step of frequency window
    freqidx : int
        Index of frequency column
    winsize : int
        Size of bp window
    winstep : int
        Step of bp window
    winidx : int
        Index of bp column
    scoreidx : int
        Index of score column
    numsnp : int
        Min # of SNPs in haplotype

    Returns
    -------
    pandas.DataFrame
        Haplotype table
    '''
    # local function
    def double_window(table,
                      leftaaf,
                      rightaaf,
                      afreqcol,
                      leftwin,
                      rightwin,
                      wincol
                     ):
        '''Grid subset of table (frequency by position)'''
        subtable=table[(table[afreqcol]>=leftaaf)&(table[afreqcol]<=rightaaf)]
        subtable=subtable[(table[wincol]>=leftwin)&(subtable[wincol]<=rightwin)]
        return subtable
    def mean(x):
        return sum(x)/len(x)
    # column names
    headernames=list(table.columns)
    wincol=headernames[winidx] # position name
    freqcol=headernames[freqidx] # frequency name
    scorecol=headernames[scoreidx] # score name
    # code
    winmin=table[wincol].min()
    winmax=table[wincol].max()
    aafmin=0
    aafmax=1
    aafrange=np.arange(aafmin,aafmax,freqstep)
    winrange=np.arange(winmin,winmax,winstep)
    windowed=dict()
    ctr=0
    for a in aafrange:
        aafleft=a
        aafright=a+freqsize
        aafmid=(aafleft+aafright)/2
        for w in winrange:
            winleft=w
            winright=w+winsize
            winmid=(winleft+winright)/2
            haplotable=double_window(table,
                                     aafleft,
                                     aafright,
                                     freqcol,
                                     winleft,
                                     winright,
                                     wincol
                                    )
            if haplotable.shape[0] >= numsnp:
                haplotable.sort_values(by=scorecol,ascending=False,inplace=True)
                score=mean(list(haplotable[scorecol])[:numsnp])
                tup=(winmid,
                     aafmid,
                     score,
                     haplotable.shape[0]
                    )
                windowed[ctr]=list(tup)
                ctr+=1
    windowed=pd.DataFrame(windowed)
    windowed=windowed.T
    windowed.columns=['POS','AAF','SCORE','SIZE']
    windowed.sort_values(by='SCORE',ascending=False,inplace=True)
    windowed.reset_index(inplace=True,drop=True)
    return windowed

# haplotype math
snptab=pd.read_csv(snpin,sep='\t')
haptab=haplotypes(snptab,
                  freqsize,
                  freqstep,
                  freqidx,
                  winsize,
                  winstep,
                  winidx,
                  scoreidx,
                  numsnp
                 )
haptab=haptab[haptab['AAF']>=low]
haptab.sort_values(by='SCORE',
                   ascending=False,
                   inplace=True,
                   ignore_index=True
                   )
haptab.to_csv(folderout+'/haplotypes.tsv',sep='\t',index=False)

# best haplotype
bestbp=haptab['POS'][0]
bestaf=haptab['AAF'][0]

f=open(folderout+'/third.best.hap.txt','w')
f.write('bp\t')
f.write(str(int(bestbp)))
f.write('\n')
f.write('frequency\t')
f.write(str(round(bestaf,4)))
f.write('\n')
f.close()

# haplotypes figure
plt.scatter(haptab['POS'],haptab['AAF'],c=haptab['SCORE'],cmap='copper_r',s=5)
plt.ylim(-0.1,1.1)
plt.colorbar(label='Haplotype z-score')
plt.ylabel('Haplotype frequency')
plt.xlabel('Position')
plt.savefig(folderout+'/third.hap.png',dpi=300)
plt.clf()

# snps figure
headernames=list(snptab.columns)
wincol=headernames[winidx] # position name
freqcol=headernames[freqidx] # frequency name
scorecol=headernames[scoreidx] # score name
plt.scatter(snptab[wincol],snptab[freqcol],c=snptab[scorecol],cmap='copper_r',s=5)
plt.ylim(-0.1,1.1)
plt.colorbar(label='SNP z-score')
plt.ylabel('SNP frequency')
plt.xlabel('Position')
plt.savefig(folderout+'/third.snp.png',dpi=300)
plt.clf()
