# Append info to region of interest
# Seth D. Temple, sdtemple@uw.edu
# May 2, 2023

# setup
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
folder,mbbuffer,mapfolder=sys.argv[1:]
filein=folder+'/excess.region.ibd.tsv'
folder1=folder+'/ibdsegs/ibdends/modified/scan/'
filepre='chr'
filesuf='.ibd.windowed.tsv.gz'
fileout=folder+'/roi.tsv'
mbbuffer=int(float(mbbuffer))*1000000

roitab=pd.read_csv(filein,sep='\t')
bpcenter=[]
cmcenter=[]

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
    # bpcenter.append(meBP)
    # cmcenter.append(meCM)
    bpcenter.append(moBP)
    cmcenter.append(moCM)

    # writing stats
    chrotit=str(int(float(rowing.CHROM)))
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

    # reading
    rowing=roitab.iloc[i]
    filein1=folder1+filepre+str(int(float(rowing.CHROM)))+filesuf
    tab=pd.read_csv(filein1, sep='\t')
    left=int(float(rowing.BPLEFT))-mbbuffer
    right=int(float(rowing.BPRIGHT))+mbbuffer
    tab=tab[(tab['BPWINDOW']>=left)&(tab['BPWINDOW']<=right)]

    # writing plots
    chrotit=str(int(float(rowing.CHROM)))
    # plotout=folder+'/plots/chr'+ chrotit + '.mean' + str(meCM) + '.cm.png'
    plotout=folder+'/plots/chr'+ chrotit + '.mode' + str(moCM) + '.cm.png'
    # plottit='Chromosome '+chrotit+', Mean '+str(meCM)
    plottit='Chromosome '+chrotit+', Mode '+str(moCM)
    plt.plot(tab['CMWINDOW'],tab['COUNT'],linewidth=3,color='k')
    plt.xlabel('centiMorgan position')
    plt.ylabel('IBD segment count')
    plt.title(plottit)
    # plt.axvline(moCM,linestyle='dotted',linewidth=2,color='tab:gray')
    # plt.axvline(mdCM,linestyle='dashed',linewidth=2,color='tab:gray')
    # plt.axvline(meCM,linestyle='solid',linewidth=2,color='tab:gray')
    plt.savefig(plotout,dpi=300)
    plt.clf()

    # writing plots
    chrotit=str(int(float(rowing.CHROM)))
    # plotout=folder+'/plots/chr'+ chrotit + '.mean' + str(meBP) + '.bp.png'
    plotout=folder+'/plots/chr'+ chrotit + '.mode' + str(moBP) + '.bp.png'
    # plottit='Chromosome '+chrotit+', Mean '+str(meBP)
    plottit='Chromosome '+chrotit+', Mode '+str(moBP)
    plt.plot(tab['BPWINDOW'],tab['COUNT'],linewidth=3,color='k')
    plt.xlabel('basepair position')
    plt.ylabel('IBD segment count')
    plt.title(plottit)
    # plt.axvline(moCM,linestyle='dotted',linewidth=2,color='tab:gray')
    # plt.axvline(mdCM,linestyle='dashed',linewidth=2,color='tab:gray')
    # plt.axvline(meBP,linestyle='solid',linewidth=2,color='tab:gray')
    plt.savefig(plotout,dpi=300)
    plt.clf()

# adding bp, cm centrality to roi table
roitab['BPCENTER']=bpcenter
roitab['CMCENTER']=cmcenter
roitab['BPLEFTCENTER']=(roitab['BPLEFT']-mbbuffer).clip(lower=1)
roitab['BPRIGHTCENTER']=roitab['BPRIGHT']+mbbuffer
# sorting, giving generic names
initcol=list(roitab.columns)
finacol=['NAME']+initcol
roitab.sort_values(by='MAXIBD',ascending=False,inplace=True)
nrow=roitab.shape[0]
roitab['NAME']=['hit'+str(i) for i in range(1,nrow+1)]
roitab=roitab[finacol]
roitab.to_csv(fileout,sep='\t',index=False)

# plotting recombination variability
J=roitab.shape[0]
for j in range(J):
    row=roitab.iloc[j]
    # chrstr=str(row.CHROM.tolist()[0])
    chrstr=str(row.CHROM)
    mapfile=mapfolder+'/chr'+chrstr+'.map'
    genmap=pd.read_csv(mapfile,sep='\t',header=None)
    # cutleft=row.BPLEFTCENTER.tolist()[0]
    cutleft=row.BPLEFTCENTER
    # cutright=row.BPRIGHTCENTER.tolist()[0]
    cutright=row.BPRIGHTCENTER
    genmapsub=genmap[(genmap[3]<cutright)&(genmap[3]>cutleft)]
    genmapsub.reset_index(inplace=True)
    cmsdiff=[]
    bpsdiff=[]
    K=genmapsub.shape[0]-1
    arr=genmapsub[2]
    arr2=genmapsub[3]
    for k in range(K):
        cmdiff=arr[k+1]-arr[k]
        bpdiff=arr2[k+1]-arr[k]
        cmsdiff.append(cmdiff)
        bpsdiff.append(bpdiff)
    diffs=np.array(cmsdiff)/np.array(bpsdiff)
    plt.scatter(genmapsub[3][:-1],diffs,s=3)
    plt.ylabel('cM / bp')
    plt.xlabel('bp position')
    # namestr=str(row.NAME.tolist()[0])
    namestr=str(row.NAME)
    chrotit=str(int(float(row.CHROM)))
    # meBP=str(int(float(row.BPCENTER)))
    # plotout=folder+'/plots/chr'+ chrotit + '.mean' + str(meBP) + '.recomb.png'
    moBP=str(int(float(row.BPCENTER)))
    plotout=folder+'/plots/chr'+ chrotit + '.mode' + str(moBP) + '.recomb.png'
    plt.title('Chromosome '+chrotit)
    plt.savefig(plotout,dpi=300)
    plt.clf()