import sys
import os
from isweep import *

# first arg: # replicates
# second arg: # bootstraps
# third arg: # diploids
# fourth arg: selection coefficient
# fifth arg: allele frequency
# sixth arg: Ne demography file
# seventh arg: output folder
# eighth arg: output file prefix

folder=sys.argv[7]
if not os.path.exists(folder):
    os.mkdir(folder)

# fixed parameter settings
# change selection coefficients if interested
lst=[1.0,2.0,3.0,4.0] # cM threshold for IBD calling

# input parameter settings
s=float(sys.argv[4])
p=float(sys.argv[5]) # allele freq
Ne=read_Ne(sys.argv[6]) # demo history
Me=sys.argv[6]

# input sim settings
nreplicates=int(float(sys.argv[1]))
nbootstraps=int(float(sys.argv[2]))
nsamples=int(float(sys.argv[3]))

# fixed sim settings
ploidy=2
msamples=ploidy*nsamples
N=msamples*(msamples-1)/2-msamples
outfile=folder+'/'+sys.argv[8]

# write to file concurrently

# f=open(outfile+'.cM.tsv','w')
f=open(outfile+'.cM.v2.tsv','w')
f.write('TRUE\tINITEST\tCORREST\tCORRLOW\tCORRUPP\tP0\tNE\tINH\tTAU\tSV\tCM\n')

sinfs=[[] for k in range(len(lst))]
for i in range(nreplicates):
    for k in range(len(lst)):
        long_ibd=lst[k]
        ab=[long_ibd,np.inf]
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
            print(lst[k],i,j) # comment out to stop stdout
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
        sinf=[s,se,sint[1],sint[0],sint[2],p,Me,0,'m',0,long_ibd]
        sinfs[k].append(sinf)
        # saving
        row=sinf
        for l in range(len(row)):
            f.write(str(row[l]))
            if l != (len(row)-1):
                f.write('\t')
            else:
                f.write('\n')

f.close()
