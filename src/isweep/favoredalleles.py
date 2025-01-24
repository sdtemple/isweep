# imports
import gzip
import numpy as np
import pandas as pd
from copy import deepcopy

### relabeling, recoding ###

def labeled_allele_frequencies(vcf, labels1):
    '''Process phased diploid VCF data to compute allele frequencies for 2 groups

    Parameters
    ----------
    vcf : str
        Name of phased diploid VCF data
    labels1 : list
        List of haplotype names in group 1

    Return
    -------
    tuple
        (basepair positions, allele frequencies group 1, allele frequencies group 0, ALT allele frequencies)
    '''

    # bring in group 1
    assert type(labels1) is list
    labhap1=labels1

    # compute allele frequencies from phased VCF
    with gzip.open(vcf, 'rt') as f:
        twochar='##'
        while twochar=='##':
            line=f.readline().strip()
            twochar=line[:2]
        header=line
        headerlist=header.split('\t')
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
    table['DELTAPRIME'] = table['DELTA'] / np.sqrt(table['AAF'] * (1 - table['AAF']))
    table.sort_values(by=['DELTAPRIME','POS'], inplace=True, ascending=False)
    table.reset_index(inplace=True,drop=True)
    if topn < table.shape[0]:
        return table.loc[:topn,]
    return table

def labeled_allele_frequencies_haploid(vcf, labels1):
    '''Process haploid VCF data to compute allele frequencies for 2 groups

    Parameters
    ----------
    vcf : str
        Name of haploid VCF data
    labels1 : list/str
        Filename for haplotype names in group 1 or a list of them

    Return
    -------
    tuple
        (basepair positions, allele frequencies group 1, allele frequencies group 0, ALT allele frequencies)
    '''

    # bring in group 1
    assert type(labels1) is list
    labhap1=labels1

    # compute allele frequencies from phased VCF
    with gzip.open(vcf, 'rt') as f:
        twochar='##'
        while twochar=='##':
            line=f.readline().strip()
            twochar=line[:2]
        header=line
        headerlist=header.split('\t')
        namdip=headerlist[9:]
        namhap1=[nam+'_1' for nam in namdip]
        hap11=[i for i in range(len(namhap1)) if namhap1[i] in labhap1]
        hap10=[i for i in range(len(namhap1)) if namhap1[i] not in labhap1]
        pos=[]
        count11=[]
        count10=[]
        itr=1
        for line in f:
            if (itr % 1000) == 0:
                print(itr)
            splitline=line.split('\t')
            infoline=splitline[:9]
            markline=splitline[9:]
            markhap1=np.array([int(float(mark[0])) for mark in markline])
            pos.append(int(float(infoline[1])))
            count11.append(markhap1[hap11].sum())
            count10.append(markhap1[hap10].sum())
            itr+=1
    denom1=len(hap11)
    denom0=len(hap10)
    count1=[(count11[i]) for i in range(len(count11))]
    count0=[(count10[i]) for i in range(len(count11))]
    freq1=np.array(count1)/denom1
    freq0=np.array(count0)/denom0
    denom=denom1+denom0
    freqm=(np.array(count1)+np.array(count0))/denom
    return pos, freq1, freq0, freqm
