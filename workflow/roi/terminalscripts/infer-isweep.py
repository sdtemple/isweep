# importing
from isweep import *
import pandas as pd
import numpy as np
import gzip
import sys

# edit these manually
inhs=['d','a','m'] # different inheritance models
svmaf=[0.05,0.02,0.01] # different standing variation
alpha=0.95
alpha1=(1-alpha)/2
alpha2=1-alpha1

# i/o
ibdct, ibdaf, fileout, nboot, cutoff, n, Ne, ploidy = sys.argv[1:]

# setting up
B=int(float(nboot))
CUTOFF=float(cutoff)
long_ibd=CUTOFF
short_ibd=long_ibd
Ne=str(Ne)
Ne=read_Ne(Ne)
ab=[CUTOFF,np.inf]
n=int(float(n))
ploidy=int(float(ploidy))
m=ploidy*n
N=m*(m-1)/2-m

# estimate allele frequency
tab=pd.read_csv(ibdaf,sep='\t')
#tab['DELTANORM']=tab['DELTA']/tab['DELTA'].sum()
p_est=(tab['AAF']*tab['DELTANORM']).sum()

# estimate selection coefficent
numTracts=0
with gzip.open(ibdct, 'rt') as f:
    for line in f:
        row=line.strip().split('\t')
        if float(row[7]) >= long_ibd:
            numTracts += 1

sinhs=[]
for inh in inhs:
    s_est = minimize_scalar(chi2_isweep,
             args=(p_est,Ne,N,(numTracts,),ab,inh),
             bounds=(0,0.5),
             method='bounded'
            ).x
    sinhs.append(s_est)

# bootstrap
sbsinhs=[[] for i in range(len(sinhs))]

for j in range(len(inhs)):
    sinh = sinhs[j]
    inh = inhs[j]
    for b in range(B):
        simdata = simulate_ibd_isweep(n,sinh,p_est,Ne,long_ibd=long_ibd,short_ibd=long_ibd,ploidy=ploidy,one_step_model=inh)
        simibd = simdata[0][0]
        sb = minimize_scalar(chi2_isweep,
             args=(p_est,Ne,N,(simibd,),ab,inh),
             bounds=(0,0.5),
             method='bounded'
            ).x
        sbsinhs[j].append(sb)

# writing in the moment
f=open(fileout,'w')
f.write('VarFreqEst\t')
f.write('SelCoefEst\t')
f.write('SelCoefLow\t')
f.write('SelCoefUpp\t')
f.write('Mendel\t')
f.write('TimeSV5Est\t')
f.write('TimeSV5Low\t')
f.write('TimeSV5Upp\t')
f.write('TimeSV2Est\t')
f.write('TimeSV2Low\t')
f.write('TimeSV2Upp\t')
f.write('TimeSV1Upp\t')
f.write('TimeSV1Est\t')
f.write('TimeSV1Low\t')
f.write('TimeDeNovoUpp\t')
f.write('TimeDeNovoEst\t')
f.write('TimeDeNovoLow\t')
f.write('TimeDeNovoUpp\n')

# correct, interval estimate
for j in range(len(inhs)):
    # writing the selection coefficient, allele frequency results
    sinh=sinhs[j]
    sbsinh=sbsinhs[j]
    inh=inhs[j]
    sl,sm,su=bootstrap_standard_bc(sinh, sbsinh, alpha1, alpha2)
    f.write(str(p_est)); f.write('\t')
    f.write(str(sm)); f.write('\t')
    f.write(str(sl)); f.write('\t')
    f.write(str(su)); f.write('\t')
    f.write(str(inh)); f.write('\t')
    svs=[[] for i in range(len(svmaf)+1)] # 1 more for de novo
    for i in range(len(svmaf)):
        # writing time until standing variation
        sv=svmaf[i]
        themeantime=when_freq(sv,sinh,p_est,Ne,one_step_model=inh,ploidy=ploidy,random_walk=False)
        for k in range(len(sbsinh)):
            # bootstrap once
            sb=sbsinh[k]
            thetime=when_freq(sv,sb,p_est,Ne,one_step_model=inh,ploidy=ploidy,random_walk=True)
            svs[i].append(thetime)
        # time estimates
        tl = np.quantile(thetime,alpha1)
        tu = np.quantile(thetime,alpha2)
        tm=themeantime
        f.write(str(tm)); f.write('\t')
        f.write(str(tl)); f.write('\t')
        f.write(str(tu)); f.write('\t')
    themeantime=when_count(1,sinh,p_est,Ne,one_step_model=inh,ploidy=ploidy,random_walk=False)
    for k in range(len(sbsinh)):
        sb=sbsinh[k]
        thetime=when_count(1,sb,p_est,Ne,one_step_model=inh,ploidy=ploidy,random_walk=True)
        svs[len(svs)].append(thetime)
    # time estimates
    tl = np.quantile(thetime,alpha1) # hall method
    tu = np.quantile(thetime,alpha2)
    tm=themeantime
    f.write(str(tm)); f.write('\t')
    f.write(str(tl)); f.write('\t')
    f.write(str(tu)); f.write('\n')
f.close()
