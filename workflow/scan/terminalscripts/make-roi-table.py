# Append info to region of interest
# Seth D. Temple, sdtemple@uw.edu
# May 2, 2023

# setup
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
folder=sys.argv[1]
filein=folder+'/excess.region.ibd.tsv'
folder1=folder+'/ibdsegs/ibdends/modified/scan/'
filepre='chr'
filesuf='.ibd.windowed.tsv.gz'
fileout=folder+'/roi.tsv'

roitab=pd.read_csv(filein,sep='\t')
bpcenter=[]
cmcenter=[]

print(roitab)

# make magic happen

# formatting
for i in range(roitab.shape[0]):

    # reading
    rowing=roitab.iloc[i]
    filein1=folder1+filepre+str(int(float(rowing.CHROM)))+filesuf
    tab=pd.read_csv(filein1, sep='\t')
    left=int(float(rowing.BPLEFT))
    right=int(float(rowing.BPRIGHT))
    tab=tab[(tab['BPWINDOW']>=left)&(tab['BPWINDOW']<=right)]
    tab['WEIGHT']=tab['COUNT']/tab['COUNT'].sum()
    tab['CUMSUM']=np.cumsum(tab['WEIGHT'])

    # central tendency
    shape0=tab.shape[0]
    mx=0
    moCM=0
    moBP=0
    print(rowing)
    for j in range(shape0):
        row=tab.iloc[j]
        if row['COUNT'] > mx:
            moCM = row['CMWINDOW']
            moBP = row['BPWINDOW']
            mx = row['COUNT']
    moBP=int(moBP)
    meCM=(tab['WEIGHT']*tab['CMWINDOW']).sum()
    meBP=(tab['WEIGHT']*tab['BPWINDOW']).sum()
    meBP=int(meBP)
    mdCM=tab[tab['CUMSUM']>=0.5]['CMWINDOW'].tolist()[0]
    mdBP=tab[tab['CUMSUM']>=0.5]['BPWINDOW'].tolist()[0]
    bpcenter.append(meBP)
    cmcenter.append(meCM)

    # writing plots
    chrotit=str(int(float(rowing.CHROM)))
    plotout=folder+'/plots/chr'+ chrotit + '.mean' + str(meBP) + '.png'
    plottit='Chromosome '+chrotit+', Mean '+str(meCM)
    plt.plot(tab['CMWINDOW'],tab['COUNT'],linewidth=3,color='k')
    plt.xlabel('centiMorgan distance')
    plt.ylabel('IBD segment count')
    plt.title(plottit)
    plt.axvline(moCM,linestyle='dotted',linewidth=2,color='tab:gray')
    plt.axvline(mdCM,linestyle='dashed',linewidth=2,color='tab:gray')
    plt.axvline(meCM,linestyle='solid',linewidth=2,color='tab:gray')
    plt.savefig(plotout,dpi=300)
    plt.clf()

    # writing stats
    statout=folder+'/stats/chr'+ chrotit + '.mean' + str(meBP) + '.tsv'
    f=open(statout,'w')
    f.write('CHROM\tMODECM\tMEDIANCM\tMEANCM\tMODEBP\tMEDIANBP\tMEANBP\n')
    f.write(chrotit); f.write('\t')
    f.write(str(moCM)); f.write('\t')
    f.write(str(mdCM)); f.write('\t')
    f.write(str(meCM)); f.write('\t')
    f.write(str(moBP)); f.write('\t')
    f.write(str(mdBP)); f.write('\t')
    f.write(str(meBP)); f.write('\n')
    f.close()

# adding bp, cm centrality to roi table
roitab['BPCENTER']=bpcenter
roitab['CMCENTER']=cmcenter
roitab.to_csv(fileout,sep='\t',index=False)
