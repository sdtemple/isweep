import sys
import os
from isweep import *

# first arg: # replicates
# second arg: # bootstraps
# third arg: # diploids
# fourth arg: allele frequency
# fifth arg: Ne demography file
# sixth arg: output file
# seventh arg: selection coefficient
# eigth arg: morgan threshold

sel=float(sys.argv[7])
lst=[sel]

# input parameter settings
p=float(sys.argv[4]) # allele freq
Ne=read_Ne(sys.argv[5]) # demo history
Me=sys.argv[5]

# input sim settings
nreplicates=int(float(sys.argv[1]))
nbootstraps=int(float(sys.argv[2]))
nsamples=int(float(sys.argv[3]))

# fixed sim settings
ploidy=2
msamples=ploidy*nsamples
long_ibd=float(sys.argv[8])
N=msamples*(msamples-1)/2-msamples
ab=[long_ibd,np.inf]
outfile=sys.argv[6]

# write to file concurrently

f=open(outfile,'w')
f.write('SELECTIONCOEFFICIENT: '); f.write(sys.argv[7]); f.write('\n')
f.write('ALLELEFREQ: '); f.write(sys.argv[4]); f.write('\n')
f.write('CMTHRESHOLD: '); f.write(sys.argv[8]); f.write('\n')
f.write('SAMPLESIZE: '); f.write(sys.argv[3]); f.write('\n')
f.write('DEMOGRAPHY: '); f.write(sys.argv[5]); f.write('\n')
f.write('NUMREPLICATES: '); f.write(sys.argv[1]); f.write('\n')
f.write('NUMBOOTSTRAPS: '); f.write(sys.argv[2]); f.write('\n')
f.write('\n') 

f.write('FORMAT\n-----\n')
f.write('ESTIMATE\tNUMIBD\n')
f.write('IBD\tIBD\t...\n')
f.write('EST\tEST\t...\n')
f.write('\n')

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
        f.write(str(se)); f.write('\t'); f.write(str(ibd)); f.write('\n')
        sboot=[]
        iboot=[]
        for j in range(nbootstraps):
            # print(lst[k],i,j) # comment out to stop stdout
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
            iboot.append(ibd)
            sboot.append(sj)
        for ib in iboot[:-1]:
            f.write(str(ib)); f.write('\t')
        f.write(str(iboot[-1])); f.write('\n')
        for sb in sboot[:-1]:
            f.write(str(sb)); f.write('\t')
        f.write(str(sboot[-1])); f.write('\n')
        f.write('\n')

f.close()
