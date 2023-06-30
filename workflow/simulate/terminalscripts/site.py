# Refine sweep site
# seth d temple, sdtemple@uw.edu
# june 30, 2023

# imports
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# inputs
(snpin, # tab-delimited snps (table 1)
 mapin, # tab-delimited w/ bp and cM (table 2)
 folderout, # folder for output
 bpidx, # column index for bp in table 1
 bpcol, # column name for bp in table 2
 cmcol, # column name for cM in table 2
 winsize, # size of sliding window
 winstep, # step of sliding window
 wincol, # column index for distance unit
 freqcol, # column index for frequency
 qrng, # sum 1... quantiles
 maxspace # max spacing filter
 )=sys.argv[1:]

bpidx=int(float(bpidx))
wincol=int(float(wincol))
freqcol=int(float(freqcol))
winsize=float(winsize)
winstep=float(winstep)
qrng=int(float(qrng))
maxspace=float(maxspace)

# functions
def q_score(table,
            afcol=1,
            dcol=-1,
            by=0.05,
            size=0.2,
            qrange=10,
            jump=0.01
           ):
    '''Compute quantile score over sliding window

    Parameters
    ----------
    table : pandas DataFrame
        SNPs with position, allele frequency
    fcol : str
        Frequency column name
    dcol : str
        Distance column name
    by : numeric
        Step for sliding window
    size : numeric
        Size for sliding window
    qrange : int
        Sum quantiles 1...
    jump : numeric
        Filter if spacing exceeds this in a window

    Returns
    -------
    pandas DataFrame
        Table for quantile scored windows
    '''
    # column names
    headernames=list(table.columns)
    pos=headernames[dcol] # position name
    aaf=headernames[afcol] # frequency name
    # spacings
    table.sort_values(by=pos,ignore_index=True,inplace=True)
    spacings=[]
    x=table[pos]
    for i in range(table.shape[0]):
        try:
            spacings.append(x[i+1]-x[i])
        except KeyError: # for last item
            spacings.append(0)
    table['SPACING']=spacings
    table['JUMP']=table['SPACING']>=jump
    # sliding windows
    winmin=table[pos].min()
    winmax=table[pos].max()-size
    winrange=np.arange(winmin,winmax,by)
    windowed=dict()
    ctr=0
    for w in winrange:
        winleft=w
        winright=w+size
        winmid=(winleft+winright)/2
        # boundaries
        subtable=table[(table[pos]>=winleft)&(table[pos]<winright)]
        # has a SNP
        sz=subtable.shape[0]
        if sz >= 1:
            # has no jump
            if subtable['JUMP'].sum() <= 0:
                sm=0
                for i in range(qrange):
                    per=(i+1)*0.01
                    val=np.quantile(subtable[aaf],per)
                    sm+=val
                sp=subtable['SPACING'].max()
                windowed[ctr]=(winmid,
                               sm,
                               sz,
                               sp
                              )
                ctr += 1
    # data frame
    windowed=pd.DataFrame(windowed)
    windowed=windowed.T
    windowed.columns=['POS',
                      'QSCORE',
                      'SIZE',
                      'MAXSPACING'
                     ]
    windowed.sort_values('QSCORE',
                         ascending=False,
                         inplace=True
                        )
    return windowed

# def centiMorgan(table,
#                 bpcol,
#                 cmcol,
#                 idx=0
#                ):
#     '''Append a column for cM distance
#
#     Parameters
#     ----------
#     table : pandas DataFrame
#         Rows are SNPs with a bp column
#     bpcol : array-like
#         bp column
#     cmcol : array-like
#         cM column (match bpcol)
#     idx : int
#         Column index with bp in table
#
#     Returns
#     -------
#     pandas DataFrame
#         Same table but with a cM column
#     '''
#     nm=list(table.columns)[idx]
#     z1=np.sort(bpcol)
#     z2=np.sort(cmcol)
#     X=[]
#     Y=table[nm]
#     for y in Y:
#         bp1=z1[z1<=y].tolist()[-1]
#         bp2=z1[z1>y].tolist()[0]
#         cm1=z2[z1<=y].tolist()[-1]
#         cm2=z2[z1>y].tolist()[0]
#         c=bp1+bp2
#         a=bp1/c
#         b=bp2/c
#         cm=cm1*a+cm2*b
#         X.append(cm)
#     table['CM']=X
#     return table

# def basepair(table,
#              bpcol,
#              cmcol,
#              idx=0
#             ):
#     '''Append a column for bp distance
#
#     Parameters
#     ----------
#     table : pandas DataFrame
#         Rows are SNPs with a bp column
#     bpcol : array-like
#         bp column
#     cmcol : array-like
#         cM column (match bpcol)
#     idx : int
#         Column index with distance in table
#
#     Returns
#     -------
#     pandas DataFrame
#         Same table but with a bp column
#     '''
#     nm=list(table.columns)[idx]
#     z1=np.sort(bpcol)
#     z2=np.sort(cmcol)
#     X=[]
#     Y=np.sort(table[nm])
#     for y in Y:
#         cm1=z2[z2<=y].tolist()[-1]
#         cm2=z2[z2>y].tolist()[0]
#         bp1=z1[z2<=y].tolist()[-1]
#         bp2=z1[z2>y].tolist()[0]
#         c=cm1+cm2
#         a=cm1/c
#         b=cm2/c
#         bp=bp1*a+bp2*b
#         X.append(bp)
#     table['BP']=X
#     return table

# compute quantile score over sliding window
snptab=pd.read_csv(snpin,sep='\t')
ibdtab=pd.read_csv(mapin,sep='\t')
# snptabcm=centiMorgan(snptab,ibdtab[bpcol],ibdtab[cmcol],bpidx)
qtab=q_score(snptab,
             freqcol,
             wincol,
             winstep,
             winsize,
             qrng,
             maxspace
             )
qtab.to_csv(folderout+'/first.tsv.gz',
            sep='\t',
            index=False
            )

# best window
best=list(qtab.iloc[0])
bestbp=best[0]

# plotting
plt.plot(qtab['POS'],
    qtab['QSCORE'],
    color='k'
    )
plt.xlabel('Score')
plt.ylabel('Position')
plt.savefig(folderout+'/first.png')

# record positions
f=open(folderout+'/first.pos.txt','w')
f.write(str(int(bestbp)))
f.write('\n')
f.close()
