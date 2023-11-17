# imports
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")

filein,fileout,windowsize,windowstep=sys.argv[1:]
windowsize=int(float(windowsize))
windowstep=int(float(windowstep))

# user-defined functions
# these will be added to isweep package

def isweep_scored_variants(table,
                           idx=(0,1,5),
                           xlabel='Position',
                           ylabel='Allele Frequency',
                           zlabel='Score',
                           cmap='cividis',
                           s=5
                          ):
    columns=list(table.columns)
    pos=columns[idx[0]]
    aaf=columns[idx[1]]
    sco=columns[idx[2]]
    x=table[pos]
    y=table[aaf]
    z=table[sco]
    plt.scatter(x,y,c=z,s=s,cmap=cmap)
    plt.colorbar(label=zlabel)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim(-0.05,1.05)
    return None

def isweep_window_scan(table,
                       column='AAFCV',
                       columnname='Coefficient of Variation'
                      ):
    x=table['POSWIN']
    y=table[column]
    plt.plot(x,y)
    plt.scatter(x,y,s=10)
    plt.ylabel(columnname)
    plt.xlabel('Position')
    return None

def isweep_windows(table,idx=(0,1),by=20_000,size=100_000):
    headernames=list(table.columns)
    pos=headernames[idx[0]]
    aaf=headernames[idx[1]]
    winmin=table[pos].min()
    winmax=table[pos].max()-size
    winrange=np.arange(winmin,winmax,by)
    windowed=dict()
    ctr=0
    for w in winrange:
        winleft=w
        winright=w+size
        winmid=(winleft+winright)/2
        subtable=table[(table[pos]>=winleft)&(table[pos]<winright)]
        if subtable.shape[0] >=1:
            stddev=subtable[aaf].std()
            avg=subtable[aaf].mean()
            perc50=np.quantile(subtable[aaf],0.5)
            perc10=np.quantile(subtable[aaf],0.1)
            perc5=np.quantile(subtable[aaf],0.05)
            perc2=np.quantile(subtable[aaf],0.02)
            perc1=np.quantile(subtable[aaf],0.01)
            cv=stddev/avg
            windowed[ctr]=(winmid,
                           avg,
                           stddev,
                           cv,
                           perc50,
                           perc10,
                           perc5,
                           perc2,
                           perc1
                          )
            ctr += 1
    windowed=pd.DataFrame(windowed)
    windowed=windowed.T
    windowed.columns=['POSWIN',
                      'AAFAVG',
                      'AAFSTDDEV',
                      'AAFCV',
                      'AAFPERC50',
                      'AAFPERC10',
                      'AAFPERC5',
                      'AAFPERC2',
                      'AAFPERC1'
                     ]
    return windowed

idx=(0,1)
tab=pd.read_csv(filein,sep='\t')
wintab=isweep_windows(tab,idx,windowstep,windowsize)
tabsort=wintab.sort_values(by='AAFCV')
p=tabsort.iloc[0]['AAFPERC1']
l=tabsort.iloc[0]['POSWIN']

# make line plot
# position by coefficient of variation
isweep_window_scan(wintab)
plt.savefig(fileout+'/isweep.window.scan.png',dpi=300)
plt.clf()

# make scatter plot for scored variants

cmap='cividis'
dot_size=5
idx=(0,1,5)
isweep_scored_variants(tab,idx,cmap=cmap,s=dot_size)
plt.savefig(fileout+'/isweep.scored.variants.png',dpi=300)
plt.clf()

# save text file
f=open(fileout+'/isweep.coef.var.txt','w')
f.write(str(round(p,5)))
f.write('\n')
f.write(str(round(l,0)))
f.write('\n')
f.close()
