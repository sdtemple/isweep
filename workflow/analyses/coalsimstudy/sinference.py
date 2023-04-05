import sys
from iSWEEP import *

# first arg: # replicates
# second arg: # bootstraps
# third arg: # diploids
# fourth arg: param fixed (s, or p0)
# fifth arg: Ne file
# sixth arg: output file name

# change
lst=[0.03,0.04,0.05] # selection coefficients
p=float(sys.argv[4]) # allele freq
Ne=read_Ne(sys.argv[5]) # demo history

nreplicates=int(float(sys.argv[1]))
nbootstraps=int(float(sys.argv[2]))
nsamples=int(float(sys.argv[3]))
ploidy=2
msamples=ploidy*nsamples
long_ibd=3.0
N=msamples*(msamples-1)/2-msamples
ab=[long_ibd,np.inf]
outfile=sys.argv[6]

sinfs=[[] for k in range(len(lst))]
for i in range(nreplicates):
    for k in range(len(lst)):
        s=lst[k]
        out=simulate_ibd_from_selective_sweep(nsamples,
                                              s,p,Ne,
                                              long_ibd,
                                              ploidy=2)
        ibd=out[0][0]
        se = minimize(chisquared_statistic_for_ibd_from_selective_sweep,
                      (0.01,),
                      args=(p,Ne,N,(ibd,),ab),
                      method='Nelder-Mead'
                     ).x[0]
        sboot=[]
        for j in range(nbootstraps):
            print(lst[k],i,j)
            out=simulate_ibd_from_selective_sweep(nsamples,
                                                  se,p,Ne,
                                                  long_ibd,
                                                  ploidy=2)
            ibd=out[0][0]
            sj = minimize(chisquared_statistic_for_ibd_from_selective_sweep,
                          (0.01,),
                          args=(p,Ne,N,(ibd,),ab),
                          method='Nelder-Mead'
                         ).x[0]
            sboot.append(sj)
        sint=bootstrap_standard_bc(se,sboot)
        sinf=[s,s,se,sint[1],sint[0],sint[2]]
        sinfs[k].append(sinf)

with open(outfile,'w') as f:
    f.write('PARAM\tTRUE\tINITEST\tCORREST\tCORRLOW\tCORRUPP\n')
    for i in range(len(sinfs)):
        sinf=sinfs[i]
        for j in range(len(sinf)):
            row=sinf[j]
            for k in range(len(row)):
                f.write(str(row[k]))
                if k != (len(row)-1):
                    f.write('\t')
                else:
                    f.write('\n')
