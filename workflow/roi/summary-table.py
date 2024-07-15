## File to summarize all of the region-specific analyses
## July 8, 2024

import pandas as pd
import sys

fileout,folder,roi,filetype,=sys.argv[1:]

tablein=pd.read_csv(folder+'/'+roi,sep='\t')
tablein=tablein[['NAME','CHROM','MAXIBD','ALPHA']]
nms = list(tablein['NAME'])
p0 = []
se = []
sl = []
su = []
ms = []
bs = []

for nm in nms:
    f=open(folder+'/'+str(nm)+'/third.best.'+filetype+'.txt')
    line=f.readline().strip().split('\t')
    bp=int(float(line[1]))
    f.close()
    bs.append(bp)
    restab=pd.read_csv(folder+'/'+str(nm)+'/results.'+filetype+'.tsv',sep='\t')
    p0.append(restab['VarFreqEst'][0])
    se.append(restab['SelCoefEst'][0])
    sl.append(restab['SelCoefLow'][0])
    su.append(restab['SelCoefUpp'][0])
    ms.append(restab['Model'][0])

tablein['LOCHAT'] = bs
tablein['PHAT'] = p0
tablein['SHAT'] = se
tablein['CONF_INT_LOW'] = sl
tablein['CONF_INT_UPP'] = su
tablein['MODEL'] = ms

tablein.to_csv(fileout,sep='\t',index=False,header=True)