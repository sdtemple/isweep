# Append info to region of interest

# setup
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


folder,cmcover,cmsmall,mbbuffer,mapfolder=sys.argv[1:]

filein=folder+'/excess.region.ibd.tsv'
folder1=folder+'/ibdsegs/ibdends/scan/'
filepre='chr'
filesuf='.ibd.windowed.tsv.gz'
fileout=folder+'/roi.tsv'
mbbuffer=float(mbbuffer)*1000000

cmcover=float(cmcover)
cmsmall=float(cmsmall)

roitab=pd.read_csv(filein,sep='\t')
bpcenter=[]
cmcenter=[]

# if there is no excess IBD make a blank file
if roitab.shape[0] == 0:
    f=open(fileout,'w')
    f.write('NAME\tCHROM\tMINIBD\tMAXIBD\tBPCENTER\tCMCENTER\tBPLEFTCENTER\tBPRIGHTCENTER\n')
    f.close()
    sys.exit()

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
    bpcenter.append(moBP)
    cmcenter.append(moCM)

    # writing stats
    chrotit=str(int(float(rowing.CHROM)))
    statout=folder+'/stats/chr'+ chrotit + '.mode' + str(moBP) + '.tsv'
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

    fig, ax = plt.subplots(2,1,figsize=(8.5,8.5))
    plt.subplots_adjust(hspace=0.5)

    # writing plots
    chrotit=str(int(float(rowing.CHROM)))
    plotout=folder+'/plots/chr'+ chrotit + '.mode' + str(moBP) + '.png'
    plottit='Chromosome '+chrotit+', Mode '+str(moCM)
    ax[0].plot(tab['CMWINDOW'],tab['COUNT'],linewidth=3,color='k')
    ax[0].set_xlabel('centiMorgan position')
    ax[0].set_ylabel('IBD segment count')
    ax[0].set_title(plottit)
    ax[0].axvline(moCM,linestyle='dotted',linewidth=2,color='tab:gray',label='Mode')
    ax[0].axvline(mdCM,linestyle='dashed',linewidth=2,color='tab:gray',label='Median')
    ax[0].axvline(meCM,linestyle='solid',linewidth=2,color='tab:gray',label='Mean')
    ax[0].legend()

    # writing plots
    chrotit=str(int(float(rowing.CHROM)))
    plottit='Chromosome '+chrotit+', Mode '+str(moBP)
    ax[1].plot(tab['BPWINDOW'],tab['COUNT'],linewidth=3,color='k')
    ax[1].set_xlabel('basepair position')
    ax[1].set_ylabel('IBD segment count')
    ax[1].set_title(plottit)
    ax[1].axvline(moBP,linestyle='dotted',linewidth=2,color='tab:gray',label='Mode')
    ax[1].axvline(mdBP,linestyle='dashed',linewidth=2,color='tab:gray',label='Median')
    ax[1].axvline(meBP,linestyle='solid',linewidth=2,color='tab:gray',label='Mode')
    fig.savefig(plotout,dpi=300)
    plt.close()
    plt.clf()

# adding bp, cm centrality to roi table
roitab['BPCENTER']=bpcenter
roitab['CMCENTER']=cmcenter

roitab=roitab[(roitab['CMRIGHT']-roitab['CMLEFT'])>=cmcover]
roitab['BPLEFTCENTER']=roitab['BPLEFT']
roitab['BPRIGHTCENTER']=roitab['BPRIGHT']
roitab.loc[(roitab['CMRIGHT']-roitab['CMLEFT'])<=cmsmall,'BPLEFTCENTER']=roitab.loc[(roitab['CMRIGHT']-roitab['CMLEFT'])<=cmsmall,'BPLEFTCENTER']-mbbuffer
roitab.loc[(roitab['CMRIGHT']-roitab['CMLEFT'])<=cmsmall,'BPRIGHTCENTER']=roitab.loc[(roitab['CMRIGHT']-roitab['CMLEFT'])<=cmsmall,'BPRIGHTCENTER']+mbbuffer
roitab['BPLEFTCENTER']=roitab['BPLEFTCENTER'].clip(lower=1)

# sorting, giving generic names
initcol=list(roitab.columns)
finacol=['NAME']+initcol
roitab.sort_values(by='MAXIBD',ascending=False,inplace=True)
nrow=roitab.shape[0]
roitab['NAME']=['hit'+str(i) for i in range(1,nrow+1)]
roitab=roitab[finacol]
roitab.to_csv(fileout,sep='\t',index=False,columns=['NAME','CHROM','MINIBD','MAXIBD','BPCENTER','CMCENTER','BPLEFTCENTER','BPRIGHTCENTER'])

# plotting recombination variability
J=roitab.shape[0]
for j in range(J):
    row=roitab.iloc[j]
    chrstr=str(row.CHROM)
    mapfile=mapfolder+'/chr'+chrstr+'.map'
    genmap=pd.read_csv(mapfile,sep='\t',header=None)
    cutleft=row.BPLEFTCENTER
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
    namestr=str(row.NAME)
    chrotit=str(int(float(row.CHROM)))
    moBP=str(int(float(row.BPCENTER)))
    plotout=folder+'/plots/chr'+ chrotit + '.mode' + str(moBP) + '.recomb.png'
    plt.title('Chromosome '+chrotit)
    plt.savefig(plotout,dpi=300)
    plt.clf()
