from iSWEEP import *
import pandas as pd
import numpy as np
import gzip

# importing
ibdct=snakemake.input.ibdlong
ibdaf=snakemake.input.ibdcomm
fileout=snakemake.output.fileout
freq_file=snakemake.input.freq

# setting up
B=int(float(snakemake.config['ISWEEPPARAMS']['NBOOT']))
CUTOFF=float(snakemake.config['ISWEEPPARAMS']['IBDCUTOFF'])
long_ibd=CUTOFF
short_ibd=float(snakemake.config['THIRDHAP']['MINOUT'])
# power=float(snakemake.config['ISWEEPPARAMS']['INVPOW'])
header=int(float(snakemake.params.header))
Ne=snakemake.config['FIXED']['iNe']
Ne=read_Ne(Ne)
ab=[CUTOFF,np.inf]
n=int(float(snakemake.config['FIXED']['SAMPSIZE']))
ploidy=int(float(snakemake.config['FIXED']['PLOIDY']))
m=ploidy*n
N=m*(m-1)/2-m

# estimate allele frequency
ktab=pd.read_csv(ibdaf,sep='\t')
p_est=(ktab.iloc[0])['AAF']
print(p_est)

# estimate selection coefficent
numTracts=0
with gzip.open(ibdct, 'rt') as f:
    for line in f:
        numTracts += 1
s_est = minimize_scalar(chi2_isweep,
         args=(p_est,Ne,N,(numTracts,),ab),
         bounds=(0,0.5),
         method='bounded'
        ).x

# bootstrap
sbs=[]
for b in range(B):
    print(b)
    simdata = simulate_ibd_isweep(n,s_est,p_est,Ne,long_ibd=long_ibd,short_ibd=short_ibd,ploidy=ploidy)
    simibd = simdata[0][0]
    sb = minimize_scalar(chi2_isweep,
         args=(p_est,Ne,N,(simibd,),ab),
         bounds=(0,0.5),
         method='bounded'
        ).x
    sbs.append(sb)

# correct, interval estimate
sl,sm,su=bootstrap_standard_bc(s_est, sbs)

# get frequency from file
f = open(freq_file, 'r')
for line in f:
    row = line.strip().split('\t')
    val = float(row[0])
p = val

# write to file
with open(fileout, 'w') as f:
    f.write('VarFreq\t')
    f.write('VarFreqEst\t')
    f.write('SelCoefEst\t')
    f.write('SelCoefLow\t')
    f.write('SelCoefMid\t')
    f.write('SelCoefUpp\n')
    # f.write('TauTime\t')
    # f.write('TauTimeEst\t')
    # f.write('TauTimeLow\t')
    # f.write('TauTimeMid\t')
    # f.write('TauTimeUpp\t')
    # f.write('AF1TimeCoef\t')
    # f.write('AF1TimeEst\t')
    # f.write('AF1TimeLow\t')
    # f.write('AF1TimeMid\t')
    # f.write('AF1TimeUpp\t')
    # f.write('AFNeaTimeCoef\t')
    # f.write('AFNeaTimeEst\t')
    # f.write('AFNeaTimeLow\t')
    # f.write('AFNeaTimeMid\t')
    # f.write('AFNeaTimeUpp\t')
    # f.write('AF5TimeCoef\t')
    # f.write('AF5TimeEst\t')
    # f.write('AF5TimeLow\t')
    # f.write('AF5TimeMid\t')
    # f.write('AF5TimeUpp\n')
    f.write(str(p))
    f.write('\t')
    f.write(str(p_est))
    f.write('\t')
    f.write(str(s_est))
    f.write('\t')
    f.write(str(sl))
    f.write('\t')
    f.write(str(sm))
    f.write('\t')
    f.write(str(su))
    f.write('\t')
    # f.write(str(new))
    # f.write('\t')
    # f.write(str(new_est))
    # f.write('\t')
    # f.write(str(cl))
    # f.write('\t')
    # f.write(str(cm))
    # f.write('\t')
    # f.write(str(cu))
    # f.write('\t')
    # f.write(str(one))
    # f.write('\t')
    # f.write(str(one_est))
    # f.write('\t')
    # f.write(str(ul))
    # f.write('\t')
    # f.write(str(um))
    # f.write('\t')
    # f.write(str(uu))
    # f.write('\t')
    # f.write(str(two))
    # f.write('\t')
    # f.write(str(two_est))
    # f.write('\t')
    # f.write(str(ml))
    # f.write('\t')
    # f.write(str(mm))
    # f.write('\t')
    # f.write(str(mu))
    # f.write('\t')
    # f.write(str(five))
    # f.write('\t')
    # f.write(str(five_est))
    # f.write('\t')
    # f.write(str(rl))
    # f.write('\t')
    # f.write(str(rm))
    # f.write('\t')
    # f.write(str(ru))
    f.write('\n')
