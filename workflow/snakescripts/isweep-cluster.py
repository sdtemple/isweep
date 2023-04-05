from isweep import *
import pandas as pd
# from copy import deepcopy
# from sklearn.cluster import KMeans

# setting up
filein=snakemake.input.filein
sklearnk=str(snakemake.config["ISWEEPPARAMS"]["SKLEARNK"])
confperc=str(snakemake.config["ISWEEPPARAMS"]["CONFPERC"])
K=int(float(sklearnk))
confperc=float(confperc)
random_state=int(float(snakemake.config["ISWEEPPARAMS"]["RANDSTATE"]))
fileout=snakemake.output.fileout
fname=snakemake.output.fignout
dpi=float(snakemake.config["ISWEEPPARAMS"]["DPI"])

# isweep clustering with kmeans
tab=pd.read_csv(filein, sep='\t')
cls=isweep_kmeans(tab,K,confperc,random_state)
plot_isweep_kmeans(tab,cls,fname,dpi)
cls.to_csv(fileout, sep='\t', index=False)

# tab=pd.read_csv(filein, sep='\t')
# tab1=deepcopy(tab)
# cc='POS'
# tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
# cc='AAF'
# tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
# cc='AAF1'
# tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
# cc='AAF0'
# tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
# cc='DELTA'
# tab1[cc]=(tab[cc]-tab[cc].mean())/tab[cc].std()
# cl=KMeans(K,random_state=random_state)
# cl.fit(tab1)
# x=tab['POS']
# y=tab['AAF']
# z=tab['DELTA']
# ys=y.std()
# ym=y.mean()
# xs=x.std()
# xm=x.mean()
# zs=z.std()
# zm=z.mean()
# cls=cl.cluster_centers_
# lbl=cl.labels_
# U=len(set(lbl))
# us=[0 for u in range(U)]
# for l in lbl:
#     us[l]+=1
# ws=[]
# for k in range(K):
#     w=[]
#     wb=cls[k][0]*xs+xm
#     wa=cls[k][1]*ys+ym
#     wd=cls[k][4]*zs+zm
#     wl=us[k]
#     ws.append((wb,wa,wd,wl))
# ws = sorted(ws, key=lambda w:w[2], reverse=True)
# haplos = pd.DataFrame(ws)
# haplos.columns=['POS','AAF','DELTA','SIZE']
