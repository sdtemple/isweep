# imports
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")

# snpfilein : str; SNP table
# ibdfilein : str; IBD table (windows)
# folderout : str; folder name of output
# windowsize : int; size of windows
# windowstep : int; size of slide

snpfilein,ibdfilein,folderout,windowsize,windowstep,quant,minaf=sys.argv[1:]
windowsize=float(windowsize)
windowstep=float(windowstep)
quant=float(quant)
minaf=float(minaf)

# user-defined functions
# these will be added to isweep package

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

# plot variants
def plot_variants(table,
                  idx=(0,1,5),
                  xlabel='Position',
                  ylabel='Allele Frequency',
                  zlabel='Score',
                  cmap='Greys',
                  s=5
                 ):
    '''Plot SNPs by position, frequency, score

    Parameters
    ------
    table : pandas.DataFrame
        Table with ranked SNPs
    idx : tuple
        Indices of basepair, frequency, score
    x,y,zlabel : str
        Name for x,y,z-axis label
    cmap : str
        Name of colormap
    s : int
        Size of scatter points

    Returns
    -------
    None
    '''
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

def plot_window(table,
                xcolumn='POS',
                xcolumnname='Position (cM)',
                ycolumn='CV',
                ycolumnname='Coefficient of Variation'
               ):
    '''Plot statistic by sliding window

    Parameters
    ------
    table : pandas.DataFrame
        Statistics table by sliding window
    x,ycolumn : str
        Statistics columns
    x,ycolumnname : str
        Plot labels

    Returns
    -------
    None
    '''
    x=table[xcolumn]
    y=table[ycolumn]
    plt.plot(x,y,color='tab:gray')
    plt.scatter(x,y,s=10,color='k')
    plt.ylabel(ycolumnname)
    plt.xlabel(xcolumnname)
    return None

def bp_window(table,
              idx=(0,1),
              by=50_000,
              size=250_000
             ):
    '''Make table with SNP variation by sliding window

    Params
    ------
    table : pandas.DataFrame
        SNP table with position, frequency
    idx : tuple
        Indices of position, frequency
    by : int
        Sliding window step (in basepair)
    size : int
        Size of window (in basepair)

    Returns
    -------
    pandas.DataFrame
        Window table with SNP variation
    '''
    headernames=list(table.columns)
    pos=headernames[idx[0]]
    aaf=headernames[idx[1]]
    winmin=table[pos].min()
    winmax=table[pos].max()-size
    winrange=np.arange(winmin,winmax,by)
    windowed=dict()
    ctr=0
    for w in winrange:
        # window bounds
        winleft=w
        winright=w+size
        winmid=(winleft+winright)/2
        subtable=table[(table[pos]>=winleft)&(table[pos]<winright)]
        if subtable.shape[0] >=1:
            # summary statistics of allele frequency in window
            stddev=subtable[aaf].std()
            avg=subtable[aaf].mean()
            cv=stddev/avg
            numsnp=subtable.shape[0]
            # record window data
            windowed[ctr]=(winmid,
                           avg,
                           stddev,
                           cv,
                           numsnp
                          )
            ctr += 1
    windowed=pd.DataFrame(windowed)
    windowed=windowed.T
    windowed.columns=['POS',
                      'AVG',
                      'STDDEV',
                      'CV',
                      'NUM'
                     ]
    return windowed

def cM_window(table,
              idx=(6,2),
              by=0.05,
              size=0.25
             ):
    '''Make table with SNP variation by sliding window

    Params
    ------
    table : pandas.DataFrame
        SNP table with position, frequency
    idx : tuple
        Indices of position, frequency
    by : int
        Sliding window step (in cM)
    size : int
        Size of window (in cM)

    Returns
    -------
    pandas.DataFrame
        Window table with SNP variation
    '''
    headernames=list(table.columns)
    pos=headernames[idx[0]]
    aaf=headernames[idx[1]]
    winmin=table[pos].min()
    winmax=table[pos].max()-size
    winrange=np.arange(winmin,winmax,by)
    windowed=dict()
    ctr=0
    for w in winrange:
        # window bounds
        winleft=w
        winright=w+size
        winmid=(winleft+winright)/2
        subtable=table[(table[pos]>=winleft)&(table[pos]<winright)]
        if subtable.shape[0] >=1:
            # summary statistics of allele frequency in window
            stddev=subtable[aaf].std()
            avg=subtable[aaf].mean()
            numsnp=subtable.shape[0]
            cv=stddev/avg
            numsnp=subtable.shape[0]
            # record window data
            windowed[ctr]=(winmid,
                           avg,
                           stddev,
                           cv,
                           numsnp
                          )
            ctr += 1
    windowed=pd.DataFrame(windowed)
    windowed=windowed.T
    windowed.columns=['POS',
                      'AVG',
                      'STDDEV',
                      'CV',
                      'NUM'
                     ]
    return windowed

def ibd_window(table1,
               table2,
               table1pos='POS',
               table2pos='CMWINDOW',
               table2ibd='COUNT'
              ):
    '''Append IBD count to window table

    Parameters
    ----------
    table1 : pandas DataFrame
        SNP variation in windows
    table2 : pandas DataFrame
        IBD count in windows
    table1pos, table2pos : str
        Column name of position
    table2ibd : str
        Column name of IBD count

    Returns
    -------
        Window table + IBD column
    '''
    X=[]
    Y=table1[table1pos]
    for y in Y:
        x=table2[(table2[table2pos]>=y)]
        x.reset_index(inplace=True,drop=True)
        X.append(x[table2ibd][0])
    table1['IBDCOUNT']=X
    return table1

def normalize(x):
    mu=x.mean()
    sig=x.std()
    return (x - mu) / sig

def score_window(table):
    '''Apply score to window table

    Score is normalized CV * IBD relative to max IBD
        CV stands for coefficient of variation

    Parameters
    ----------
    table : pandas.DataFrame
        Window table with SNP variation + IBD column

    Returns
    -------
    pandas.DataFrame
        Table with a score column
    '''
    relibd=table['IBDCOUNT']/table['IBDCOUNT'].max()
    table['SCORE']=normalize(-table['CV'])*relibd
    return table

def subset_window(table,
                  center,
                  size=0.25,
                  pos='CM'
                 ):
    '''Subset SNPs into a window

    Parameter
    ---------
    table : pandas DataFrame
        SNP table
    center : float
        In bp or cM
    size : numeric
        In bp or cM
    pos : str
        Column name
    '''
    pm=size/2
    left=center-pm
    right=center+pm
    subtable=table[(table[pos]>=left)&(table[pos]<=right)]
    return subtable

def isweep_frequency(table,
                     quant=0.25,
                     minaf=0.15,
                     score='ZSCORE',
                     aaf='AAF'
                    ):
    '''Estimate allele frequency

    Parameters
    ----------
    table : pandas.DataFrame
        SNP table with scores, frequency
    quant : float
        SNPs less than this quantile used
    minaf : float
        SNPs greater than this frequency used
    score : str
        Score column name
    aaf : str
        Allele frequency column name

    Returns
    -------
    float
        An allele frequency
    '''
    subtable=table[table[aaf]>=minaf]
    cut=np.quantile(table[aaf],quant)
    aggtable=subtable[subtable[aaf]<=cut]
    aggtable['WEIGHT']=aggtable[score]/aggtable[score].sum()
    estimate=(aggtable[aaf]*aggtable['WEIGHT']).sum()
    return estimate

idx=(6,2)
snptab=pd.read_csv(snpfilein,sep='\t')
ibdtab=pd.read_csv(ibdfilein,sep='\t')
snptabcm=centiMorgan(snptab,ibdtab['BPWINDOW'],ibdtab['CMWINDOW'],0)
wintab=cM_window(snptabcm,idx,windowstep,windowsize)
fintab=ibd_window(wintab,ibdtab)
scotab=score_window(fintab)
try:
    # line plot : position by coefficient of variation
    plot_window(scotab,ycolumn='CV',ycolumnname='Coefficient of Variation')
    plt.savefig(folderout+'/isweep.window.cv.png',dpi=300)
    plt.clf()
    # line plot : position by score
    plot_window(scotab,ycolumn='SCORE',ycolumnname='Score')
    plt.savefig(folderout+'/isweep.window.score.png',dpi=300)
    plt.clf()
except:
    pass
# save unsorted
scotab.to_csv(folderout+'/isweep.window.score.unsorted.tsv',sep='\t')
# sort
scotab.sort_values(by='SCORE',inplace=True,ascending=False)
scotab.to_csv(folderout+'/isweep.window.score.sorted.tsv',sep='\t')
l=scotab.iloc[0]['POS']
b=ibdtab[ibdtab['CMWINDOW']>=l]['BPWINDOW'].tolist()[0]

try:
    # make scatter plot for scored variants
    cmap='Greys'
    dot_size=5
    idx=(6,1,5)
    plot_variants(snptabcm,idx,cmap=cmap,s=dot_size)
    plt.savefig(folderout+'/isweep.variants.png',dpi=300)
    plt.clf()
except:
    pass

try:
    # ibd count line plot
    x=ibdtab['CMWINDOW']
    y=ibdtab['COUNT']
    plt.plot(x,y,color='black')
    plt.ylabel('IBD Segment Count')
    plt.xlabel('Position (cM)')
    plt.savefig(folderout+'/isweep.ibd.png',dpi=300)
    plt.clf()
except:
    pass

# save locus text file
f=open(folderout+'/isweep.locus.txt','w')
f.write(str(round(b,0)))
f.write('\t')
f.write(str(round(l,4)))
f.write('\n')
f.close()

# save frequency estimate
subtable=subset_window(snptabcm,l,windowstep,'CM')
subtable['ZSCORE']=(subtable['AAF1']-subtable['AAF0'])/np.sqrt(subtable['AAF']*(1-subtable['AAF']))
p=isweep_frequency(subtable,quant,minaf)
f=open(folderout+'/isweep.freq.txt','w')
f.write(str(round(p,4)))
f.write('\n')
f.close()
# this is lazy --->
# could re-run cluster procedure about new locus
