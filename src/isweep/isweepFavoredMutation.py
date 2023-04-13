import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from copy import deepcopy

def labeled_allele_frequencies(vcf, labels1):
    '''Process phased diploid VCF data to compute allele frequencies for 2 groups

    Parameters
    ----------
    vcf : str
        Name of phased diploid VCF data
    labels1 : list/str
        Filename for haplotype names in group 1 or a list of them

    Return
    -------
    tuple
        (basepair positions, allele frequencies group 1, allele frequencies group 0, ALT allele frequencies)
    '''

    # bring in group 1
    if type(labels1) is list:
        labhap1=labels1
    elif type(labels1) is str:
        labhap1=[]
        with gzip.open(label1file, 'rt') as g:
            g.readline()
            for line in g:
                splitline=line.split('\t')
                labhap1.append(splitline[0])
    else:
        raise 'labels should be lists or file names'

    # compute allele frequencies from phased VCF
    with gzip.open(vcf, 'rt') as f:
        twochar='##'
        while twochar=='##':
            line=f.readline().strip()
            twochar=line[:2]
        header=line
        headerlist=header.split('\t')
        numdip=len(header)-9
        namdip=headerlist[9:]
        namhap1=[nam+'_1' for nam in namdip]
        namhap2=[nam+'_2' for nam in namdip]
        hap11=[i for i in range(len(namhap1)) if namhap1[i] in labhap1]
        hap21=[i for i in range(len(namhap2)) if namhap2[i] in labhap1]
        hap10=[i for i in range(len(namhap1)) if namhap1[i] not in labhap1]
        hap20=[i for i in range(len(namhap2)) if namhap2[i] not in labhap1]
        pos=[]
        count11=[]
        count21=[]
        count10=[]
        count20=[]
        itr=1
        for line in f:
            if (itr % 1000) == 0:
                print(itr)
            splitline=line.split('\t')
            infoline=splitline[:9]
            markline=splitline[9:]
            markhap1=np.array([int(float(mark[0])) for mark in markline])
            markhap2=np.array([int(float(mark[2])) for mark in markline])
            pos.append(int(float(infoline[1])))
            count11.append(markhap1[hap11].sum())
            count21.append(markhap2[hap21].sum())
            count10.append(markhap1[hap10].sum())
            count20.append(markhap2[hap20].sum())
            itr+=1
    denom1=len(hap11)+len(hap21)
    denom0=len(hap10)+len(hap20)
    count1=[(count11[i] + count21[i]) for i in range(len(count11))]
    count0=[(count10[i] + count20[i]) for i in range(len(count11))]
    freq1=np.array(count1)/denom1
    freq0=np.array(count0)/denom0
    denom=denom1+denom0
    freqm=(np.array(count1)+np.array(count0))/denom
    return pos, freq1, freq0, freqm

def putative_allele_frequencies(freq1, freq0, freqm):
    '''Recode allele frequencies to be putative adaptive allele frequencies

    Parameters
    ----------
    freq1 : array-like
    freq0 : array-like
    freqm : array-like
        ALT allele frequencies for group 1
        ALT allele frequencies for group 0
        ALT allele frequencies for both groups

    Return
    ------
    tuple
        array-like putative adaptive allele frequencies
        where putative adaptive is the allele more frequent in group 1 than 0
    '''
    u = []
    v = []
    w = []
    for i in range(len(freq1)):
        y = freq1[i]
        x = freq0[i]
        z = freqm[i]
        if y - z < 0:
            x = 1-x
            y = 1-y
            z = 1-z
        u.append(x)
        v.append(y)
        w.append(z)
    return np.array(v), np.array(u), np.array(w)

def format_allele_table(pos, freq1, freq0, freqm, topn=np.inf):
    '''Create a pandas DataFrame from grouped allele frequencies

    Parameters
    ----------
    pos : array-like
    freq1 : array-like
    freq0 : array-like
    freqm : array-like
        Allele frequencies for group 1, group 0, and both groups
    topn : int
        Number of top positions to study

    Return
    ------
    pandas DataFrame
        Table of adaptive allele frequencies by group
    '''
    v,u,w=putative_allele_frequencies(freq1,freq0,freqm)
    table = {'POS':pos, 'AAF':w, 'AAF1':v, 'AAF0':u, 'DELTA':v-u}
    table = pd.DataFrame(table)
    table.sort_values(by=['DELTA','POS'], inplace=True, ascending=False)
    table.reset_index(inplace=True,drop=True)
    if topn < table.shape[0]:
        return table.loc[:topn,]
    return table

# old

# def plot_Delta_by_position(table, vline=np.inf, colors=['tab:purple','tab:blue','tab:cyan','0.8'], tops=[10,20], dotsize=15, title='', loc='upper right'):
#     '''Plot delta allele frequencies to study favored mutation in iSWEEP inference
#
#     Parameters
#     ----------
#     table : pandas DataFrame
#         VCF data sorted by Delta derived allele frequency
#     vline : float
#         A basepair position for reference
#     colors : list
#         Must be 4 colors (default is ['tab:purple','tab:blue','tab:cyan','0.8'])
#     tops : array-like
#         First and second color groupings, i.e. Top 1, 1-10, and 11-20 is [10,20]
#     dotsize : int
#         Default is 15
#     title : str
#     loc : str
#         Positioning of legend
#
#     Return
#     ------
#     None
#         Makes a summary plot
#     '''
#     # try shades
#     # colors=['0','0.5','0.75','0.9']
#     a=tops[0]
#     b=tops[1]
#     labels=['Top 1','1-'+str(a),str(a)+'-'+str(b),str(b)+'+']
#     assert table.shape[0] >= 10
#     assert len(colors) == 4
#     assert len(tops) == 2
#     assert tops[0] <= tops[1]
#     xstr='POS'
#     ystr='DELTA'
#     table.sort_values(['DELTA'],inplace=True,ascending=False)
#     table.reset_index(inplace=True,drop=True)
#     subtab=table.loc[b:,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[3], label=labels[3])
#     subtab=table.loc[a:b,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[2], label=labels[2])
#     subtab=table.loc[1:a,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[1], label=labels[1])
#     subtab=table.loc[0,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[0], label=labels[0])
#     plt.axvline(x=vline,color='tab:red',ls='dotted')
#     plt.ylim(-0.1,1.1)
#     plt.xlim(table['POS'].min()-250000,table['POS'].max()+250000)
#     plt.xlabel('Relative position')
#     plt.ylabel('Delta allele frequency')
#     plt.title(title)
#     plt.legend(loc=loc)
#     return None

# def plot_Delta_by_frequency(table, vline=np.inf, colors=['tab:purple','tab:blue','tab:cyan','0.8'], tops=[10,20], dotsize=15,title='', loc='upper right'):
#     '''Plot Delta by adaptive allele frequency to study favored mutation in iSWEEP inference
#
#     Parameters
#     ----------
#     table : pandas DataFrame
#         VCF data sorted by Delta derived allele frequency
#     vline : float
#         An allele frequency for reference
#     colors : list
#         Must be 3 colors (default is ['tab:purple','tab:blue','tab:cyan','0.8'])
#     tops : array-like
#         First and second color groupings, i.e. Top 10 and 11 - 20 is [10,20]
#     dotsize : int
#         Default is 15
#     title : str
#     loc : str
#         Positioning of legend
#
#     Return
#     ------
#     None
#         Makes a summary plot
#     '''
#     # try shades
#     # colors=['0','0.5','0.75','0.9']
#     a=tops[0]
#     b=tops[1]
#     labels=['Top 1','1-'+str(a),str(a)+'-'+str(b),str(b)+'+']
#     assert table.shape[0] >= 10
#     assert len(colors) == 4
#     assert len(tops) == 2
#     assert tops[0] <= tops[1]
#     xstr='AAF'
#     ystr='DELTA'
#     table.sort_values(['DELTA'],inplace=True,ascending=False)
#     table.reset_index(inplace=True,drop=True)
#     subtab=table.loc[b:,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[3], label=labels[3])
#     subtab=table.loc[a:b,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[2], label=labels[2])
#     subtab=table.loc[1:a,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[1], label=labels[1])
#     subtab=table.loc[0,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[0], label=labels[0])
#     plt.axvline(x=vline,color='tab:red',ls='dotted')
#     plt.ylim(-0.1,1.1)
#     plt.xlim(-0.1,1.1)
#     plt.xlabel('Adaptive allele frequency')
#     plt.ylabel('Delta allele frequency')
#     plt.title(title)
#     plt.legend(loc=loc)
#     return None

# def plot_frequency_by_position(table, vline=np.inf, colors=['tab:purple','tab:blue','tab:cyan','0.8'], tops=[10,20], dotsize=15,title='', loc='upper right', method='isweep'):
#     '''Plot adaptive allele frequency by position to study favored mutation in iSWEEP inference
#
#     Parameters
#     ----------
#     table : pandas DataFrame
#         VCF data sorted by Delta derived allele frequency
#     vline : float
#         An allele frequency for reference
#     colors : list
#         Must be 3 colors (default is ['tab:purple','tab:blue','tab:cyan','0.8'])
#     tops : array-like
#         First and second color groupings, i.e. Top 10 and 11 - 20 is [10,20]
#     dotsize : int
#         Default is 15
#     title : str
#     loc : str
#         Positioning of legend
#     method : str
#         'isweep' or 'safe'
#
#     Return
#     ------
#     None
#         Makes a summary plot
#     '''
#     # try shades
#     # colors=['0','0.5','0.75,'0.9']
#     a=tops[0]
#     b=tops[1]
#     labels=['Top 1','1-'+str(a),str(a)+'-'+str(b),str(b)+'+']
#     assert table.shape[0] >= 10
#     assert len(colors) == 4
#     assert len(tops) == 2
#     assert tops[0] <= tops[1]
#     if method=='isweep':
#         xstr='POS'
#         ystr='AAF'
#         table.sort_values(['DELTA'],inplace=True,ascending=False)
#         table.reset_index(inplace=True,drop=True)
#     elif method=='isafe':
#         xstr='POS'
#         ystr='DAF'
#         table.sort_values(['iSAFE'],inplace=True,ascending=False)
#         table.reset_index(inplace=True,drop=True)
#     else:
#         raise 'Invalid method string'
#     subtab=table.loc[b:,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[3], label=labels[3])
#     subtab=table.loc[a:b,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[2], label=labels[2])
#     subtab=table.loc[1:a,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[1], label=labels[1])
#     subtab=table.loc[0,]
#     plt.scatter(subtab[xstr],subtab[ystr], s=dotsize, color=colors[0], label=labels[0])
#     plt.axvline(x=vline,color='tab:red',ls='dotted')
#     plt.ylim(-0.1,1.1)
#     plt.xlim(table[xstr].min()-250000,table[xstr].max()+250000)
#     plt.xlabel('Relative position')
#     plt.ylabel('Adaptive allele frequency')
#     plt.title(title)
#     plt.legend(loc=loc)
#     return None

# def inverse_rank_weighted_mean_std(table, topn=np.inf, power=1, method='isweep'):
#     '''Compute inverse rank weighted mean, standard deviation for adaptive allele frequency
#
#     Parameters
#     ----------
#     table : pandas DataFrame
#         VCF data sorted by delta derived allele frequency
#     topn : int
#         Number of top positions to study
#     power : float
#         Scalar to serve as exponent for inverse rank
#     method : str
#         'isweep' or 'safe'
#
#     Returns
#     -------
#     tuple
#         Mean, standard deviation estimates
#     '''
#     if method=='isweep':
#         ystr='AAF'
#         table.sort_values(['DELTA'],inplace=True,ascending=False)
#         table.reset_index(inplace=True,drop=True)
#     elif method=='isafe':
#         ystr='DAF'
#         table.sort_values(['iSAFE'],inplace=True,ascending=False)
#         table.reset_index(inplace=True,drop=True)
#     ranks = np.array([i+1 for i in range(table.shape[0])])
#     topn = min(table.shape[0], topn)
#     toprk = ranks[:topn]
#     invrk = 1 / toprk
#     scark = np.power(invrk, power)
#     norrk = scark / scark.sum()
#     freqs = table[ystr][:topn]
#     estim = (freqs * norrk).sum()
#     mome2 = (np.power(freqs,2) * norrk).sum()
#     return estim, (mome2 - estim ** 2) ** 0.5

# def violinplot_rankers(results, method_names):
#     '''Make violin plot to compare rank methods
#
#     Parameters
#     ----------
#     results : tuple
#         Array with ranks for each replicate study
#     method_names : list
#         Names of rank methods
#
#     Return
#     ------
#     None
#         Create a summary plot
#     '''
#     assert len(results) == len(method_names)
#     J = len(results[0])
#     K = len(method_names)
#     methods = []
#     for k in range(K):
#         method = [method_names[k] for j in range(J)]
#         methods += method
#     together=np.concatenate(results)
#     ax=sns.violinplot(methods, together)
#     ax.set_ylabel('Rank of true mutation')
#     return None

# def isweep_kmeans(tab,K=10,conf=0.95,random_state=None):
#     '''Perform k-means clustering to group by IBD selection signal
#
#     Parameters
#     ----------
#     tab : pandas DataFrame
#         Summary table of ranked mutations
#     K : int
#         # clusters
#     conf : float
#         like 0.XX confidence interval -> sets quantiles
#     random_state : int
#         See sklearn.cluster.KMeans
#
#     Returns
#     -------
#     pandas DataFrame
#         Summary table of ranked groups
#     '''
#     # standardizing / normalizing
#     tab1=deepcopy(tab)
#     cc='POS'
#     tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
#     cc='AAF'
#     tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
#     cc='AAF1'
#     tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
#     cc='AAF0'
#     tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
#     cc='DELTA'
#     tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
#     # k-means clustering
#     cl=KMeans(K,random_state=random_state)
#     cl.fit(tab1)
#     # setting up
#     cls=cl.cluster_centers_
#     lbl=cl.labels_
#     X=len(set(lbl))
#     vs=[0 for x in range(X)]
#     xs=[[] for x in range(X)]
#     ys=[[] for x in range(X)]
#     zs=[[] for x in range(X)]
#     alpha=1-conf
#     alpha1=alpha/2
#     alpha2=1-alpha1
#     # gathering scatter for each cluster
#     for i in range(len(lbl)):
#         l=lbl[i]
#         row=tab.iloc[i]
#         vs[l]+=1
#         ys[l].append(row['AAF'])
#         xs[l].append(row['POS'])
#         zs[l].append(row['DELTA'])
#     # summarizing
#     ws=[]
#     for k in range(K):
#         w=[]
#         # medians
#         wx=np.quantile(np.array(xs[k]),0.5)
#         wy=np.quantile(np.array(ys[k]),0.5)
#         wz=np.quantile(np.array(zs[k]),0.5)
#         wv=vs[k]
#         # lower quantile
#         wx1=np.quantile(np.array(xs[k]),alpha1)
#         wy1=np.quantile(np.array(ys[k]),alpha1)
#         wz1=np.quantile(np.array(zs[k]),alpha1)
#         # upper quantile
#         wx2=np.quantile(np.array(xs[k]),alpha2)
#         wy2=np.quantile(np.array(ys[k]),alpha2)
#         wz2=np.quantile(np.array(zs[k]),alpha2)
#         # to list
#         ws.append((wx,wy,wz,wv,wx1,wx2,wy1,wy2,wz1,wz2))
#     # sorting by delta statistic
#     ws = sorted(ws, key=lambda w:w[2], reverse=True)
#     haplos = pd.DataFrame(ws)
#     haplos.columns=['POS','AAF','DELTA','SIZE',
#                     'POSLOW','POSUPP', # position
#                     'AAFLOW','AAFUPP', # allele frequency
#                     'DELTALOW','DELTAUPP' # rank diff statistic
#                    ]
#     return haplos
#
# def plot_isweep_kmeans(tab,cls,fname='example.png',dpi=300):
#     '''Plot results of summary ranking and clustering
#
#     Parameters
#     ----------
#     tab : pandas DataFrame
#         Summary table of ranked mutations
#     cls : pandas DataFrame
#         Summary table of ranked groups
#     fname : str
#     dpi : float
#         See matplotlib.pyplot.savefig
#
#     Returns
#     -------
#     None
#         Makes a plot
#     '''
#     # main plotting
#     sns.scatterplot(x = "POS",
#                     y = "AAF",
#                     data = tab,
#                     hue = "DELTA",
#                     palette = "icefire",
#                     size='DELTA',
#                     sizes=(50,100))
#     # axes, labels, legend
#     plt.ylim(-0.2,1.2)
#     plt.legend(loc='upper center',
#                bbox_to_anchor=(0.5, 1.05),
#                ncol=10,
#                fancybox=True,
#                shadow=True,
#                title='Delta')
#     plt.xlabel('Position')
#     plt.ylabel('Allele frequency')
#     # overplotting
#     for i in range(cls.shape[0]):
#         row=cls.iloc[i]
#         plt.scatter(row[0],row[1],c='tab:purple',edgecolor='black',s=100)
#     # saving
#     plt.savefig(fname=fname,dpi=dpi)
#     return None
