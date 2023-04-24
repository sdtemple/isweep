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
# change allele frequency if interested
# lst=[-0.1,
#      -0.05,
#      0,
#      0.05,
#      0.1
#      ] # p shifts
lst=[-0.2,-0.02,0.02,0.2] # p shift

# input parameter settings
s=float(sys.argv[4]) # selection coefficient
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
long_ibd=3.0
N=msamples*(msamples-1)/2-msamples
ab=[long_ibd,np.inf]

outfile=folder+'/'+sys.argv[8]

# write to file concurrently

f=open(outfile+'.misp0.tsv','w')
f.write('TRUE\tINITEST\tCORREST\tCORRLOW\tCORRUPP\tP0\tNE\tINH\tTAU\tSV\tCM\n')

sinfs=[[] for k in range(len(lst))]
for i in range(nreplicates):
    for k in range(len(lst)):
        pk=p-lst[k]
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
            args=(pk,Ne,N,(ibd,),ab),
            bounds=(0,0.5),
            method='bounded'
        ).x
        sboot=[]
        for j in range(nbootstraps):
            print(lst[k],i,j) # comment out to stop stdout
            out=simulate_ibd_isweep(
                nsamples,
                se,
                pk,
                Ne,
                long_ibd,
                long_ibd,
                ploidy=2
            )
            ibd=out[0][0]
            sj = minimize_scalar(
                chi2_isweep,
                args=(pk,Ne,N,(ibd,),ab),
                bounds=(0,0.5),
                method='bounded'
            ).x
            sboot.append(sj)
        sint=bootstrap_standard_bc(se,sboot)
        # sinf=[s,se,sint[1],sint[0],sint[2],round(pk,3),Me,0,'m',0,long_ibd]
        sinf=[s,se,sint[1],sint[0],sint[2],round(pk,3),Me,'m',0,0,long_ibd]
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
