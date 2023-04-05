import sys
import os
from isweep import *

# first arg: # replicates
# second arg: # bootstraps
# third arg: # diploids
# fourth arg: allele frequency
# fifth arg: Ne demography file
# sixth arg: output folder
# seventh arg: output file prefix

folder=sys.argv[6]
if not os.path.exists(folder):
    os.mkdir(folder)

# fixed parameter settings
# change selection coefficients if interested
lst=[0.02,0.03,0.04,0.05,0.06] # selection coefficients

# input parameter settings
p=float(sys.argv[4]) # allele freq
Ne=read_Ne(sys.argv[5]) # demo history

# input sim settings
nreplicates=int(float(sys.argv[1]))
nbootstraps=int(float(sys.argv[2]))
nsamples=int(float(sys.argv[3]))

# fixed sim settings
ploidy=2
msamples=ploidy*nsamples
long_ibd=3.0
N=msamples*(msamples-1)/2-msamples
ab=[long_ibd,np.inf]
outfile=folder+'/'+sys.argv[7]

sinfs=[[] for k in range(len(lst))]
for i in range(nreplicates):
    for k in range(len(lst)):
        s=lst[k]
        out=simulate_ibd_isweep(
            nsamples,
            s,
            p,
            Ne,
            long_ibd,
            long_ibd,
            ploidy=2
        )
        ibd=out[0][0]
        se = minimize_scalar(
            chi2_isweep,
            args=(p,Ne,N,(ibd,),ab),
            bounds=(0,0.5),
            method='bounded'
        ).x
        sboot=[]
        for j in range(nbootstraps):
            print(lst[k],i,j)
            out=simulate_ibd_isweep(
                nsamples,
                se,
                p,
                Ne,
                long_ibd,
                long_ibd,
                ploidy=2
            )
            ibd=out[0][0]
            sj = minimize_scalar(
                chi2_isweep,
                args=(p,Ne,N,(ibd,),ab),
                bounds=(0,0.5),
                method='bounded'
            ).x
            sboot.append(sj)
        sint=bootstrap_standard_bc(se,sboot)
        sinf=[s,s,se,sint[1],sint[0],sint[2]]
        sinfs[k].append(sinf)

with open(outfile+'.tsv','w') as f:
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
