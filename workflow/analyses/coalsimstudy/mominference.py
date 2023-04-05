import sys
from iSWEEP import *

# first arg: # replicates
# second arg: # bootstraps
# third arg: # diploids
# fourth arg: param fixed (s, p, or Ne)
# fifth arg: param fixed (s, p, or Ne)
# sixth arg: output file prefix
# seventh arg: example bins file

# change
lst=[0.02,0.03,0.04,0.05,0.06] # selection coefficients
p=float(sys.argv[4]) # allele freq
Ne=read_Ne(sys.argv[5]) # demo history
bn=read_bins(sys.argv[7])

nreplicates=int(float(sys.argv[1]))
nbootstraps=int(float(sys.argv[2]))
nsamples=int(float(sys.argv[3]))
ploidy=2
msamples=ploidy*nsamples
long_ibd=3.0
N=msamples*(msamples-1)/2-msamples
outfile=sys.argv[6]

# goodness of fit to length distribution
bn.append(np.inf)
ab=bn
sinfs=[[] for k in range(len(lst))]
for i in range(nreplicates):
    for k in range(len(lst)):
        s=lst[k]
        out=simulate_ibd_from_selective_sweep(nsamples,
                                              s,p,Ne,
                                              long_ibd,
                                              ploidy=ploidy)
        ibd=big_format_distribution(out[0][2],out[0][4])
        ibd=bin_ibd_segments(ibd,ab)
        print(N)
        se = minimize(chisquared_statistic_for_ibd_from_selective_sweep, 
                      (0.01,), 
                      args=(p,Ne,N,ibd,ab), 
                      bounds=[(0,0.5)], 
                      method='Nelder-Mead'
                     ).x[0]
        sboot=[]
        for j in range(nbootstraps):
            print(lst[k],i,j)
            out=simulate_ibd_from_selective_sweep(nsamples,
                                                  se,p,Ne,
                                                  long_ibd,
                                                  ploidy=ploidy)
            ibd=big_format_distribution(out[0][2],out[0][4])
            ibd=np.array(bin_ibd_segments(ibd,ab))
            print(N)
            sj = minimize(chisquared_statistic_for_ibd_from_selective_sweep,
                          (0.01,),
                          args=(p,Ne,N,ibd,ab),
                          method='Nelder-Mead',
                          bounds=[(0,0.5)]
                         ).x[0]
            sboot.append(sj)
        sint=bootstrap_standard_bc(se,sboot)
        sinf=['GOOD',s,se,sint[1],sint[0],sint[2]]
        sinfs[k].append(sinf)

with open(outfile+'.good.tsv','w') as f:
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

# method of moments to ibd count
ab=[long_ibd,np.inf]
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
                      method='Nelder-Mead',
                      bounds=[(0,0.5)]
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
                          method='Nelder-Mead',
                          bounds=[(0,0.5)]
                         ).x[0]
            sboot.append(sj)
        sint=bootstrap_standard_bc(se,sboot)
        sinf=['MOME',s,se,sint[1],sint[0],sint[2]]
        sinfs[k].append(sinf)

with open(outfile+'.mome.tsv','w') as f:
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
