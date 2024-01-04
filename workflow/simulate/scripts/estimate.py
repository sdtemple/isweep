# Purpose: estimate s from ibd data
# Input: ibd data, p_est, cutoff, n, Ne, ploidy, nboot
# Output: s_est, s_low, s_upp

# importing
from isweep import *
import pandas as pd
import numpy as np
import gzip
import sys

# i/o
fileout, ibdct, p_est, nboot, cutoff, n, Ne, ploidy = sys.argv[1:]

# edit these manually
inh='a'
alpha=0.05
alpha2=(1-alpha/2)
alpha1=alpha/2
svmaf=[0.05,0.02] # different standing variation

# setting up
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
p_est=float(p_est)
numTracts=int(float(ibdct))
nboot=int(float(nboot))

s_est = minimize_scalar(chi2_isweep,
            args=(p_est,Ne,N,(numTracts,),ab,inh),
            bounds=(0,0.5),
            method='bounded'
        ).x

# bootstrap
sbs = []

for b in range(nboot):
    print(b)
    simdata = simulate_ibd_isweep(n,
                                  s_est,
                                  p_est,
                                  Ne,
                                  long_ibd=long_ibd,
                                  short_ibd=long_ibd,
                                  ploidy=ploidy,
                                  one_step_model=inh
                                  )
    simibd = simdata[0][0]
    sb = minimize_scalar(chi2_isweep,
            args=(p_est,Ne,N,(simibd,),ab,inh),
            bounds=(0,0.5),
            method='bounded'
        ).x
    sbs.append(sb)

with open(fileout,'w') as f:

    # writing in the moment
    f.write('VarFreqEst\t')
    f.write('SelCoefEst\t')
    f.write('SelCoefLow\t')
    f.write('SelCoefUpp\t')
    f.write('Model\t')
    f.write('TimeSV5Est\t')
    f.write('TimeSV5Low\t')
    f.write('TimeSV5Upp\t')
    f.write('TimeSV2Est\t')
    f.write('TimeSV2Low\t')
    f.write('TimeSV2Upp\t')
    f.write('TimeDeNovoEst\t')
    f.write('TimeDeNovoLow\t')
    f.write('TimeDeNovoUpp\n')

    # correct, interval estimate
    # estimator is unbiased from beginning
    # sensitivity near zero boundary, small s <= 0.02
    # sl,sm,su=bootstrap_standard_bc(sinh, sbsinh, alpha1, alpha2)
    sl,sm,su=bootstrap_standard(s_est, sbs, alpha1, alpha2)
    f.write(str(p_est)); f.write('\t')
    f.write(str(sm)); f.write('\t')
    f.write(str(sl)); f.write('\t')
    f.write(str(su)); f.write('\t')
    f.write(str(inh)); f.write('\t')

    svs=[[] for i in range(len(svmaf)+1)] # 1 more for de novo
    for i in range(len(svmaf)):
        # writing time until standing variation
        sv=svmaf[i]
        themeantime=when_freq(sv,
                              s_est,
                              p_est,
                              Ne,
                              one_step_model=inh,
                              ploidy=ploidy,
                              random_walk=False
                              )
        for k in range(len(sbs)):
            # bootstrap once
            sb=sbs[k]
            thetime=when_freq(sv,
                              sb,
                              p_est,
                              Ne,
                              one_step_model=inh,
                              ploidy=ploidy,
                              random_walk=True
                              )
            svs[i].append(thetime)
        # time estimates
        tl = np.quantile(svs[i],alpha1)
        tu = np.quantile(svs[i],alpha2)
        # tm=themeantime
        tm = np.quantile(svs[i],0.5)
        f.write(str(tm)); f.write('\t')
        f.write(str(tl)); f.write('\t')
        f.write(str(tu)); f.write('\t')
    
    themeantime=when_count(1,
                           s_est,
                           p_est,
                           Ne,
                           one_step_model=inh,
                           ploidy=ploidy,
                           random_walk=False
                           )
    for k in range(len(sbs)):
        sb=sbs[k]
        thetime=when_count(1,
                           sb,
                           p_est,
                           Ne,
                           one_step_model=inh,
                           ploidy=ploidy,
                           random_walk=True
                           )
        svs[-1].append(thetime)
    # time estimates
    tl = np.quantile(svs[-1],alpha1) # hall method
    tu = np.quantile(svs[-1],alpha2)
    # tm=themeantime
    tm = np.quantile(svs[-1],0.5)
    f.write(str(tm)); f.write('\t')
    f.write(str(tl)); f.write('\t')
    f.write(str(tu)); f.write('\n')
