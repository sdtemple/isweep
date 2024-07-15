import sys
import pandas as pd
filein,fileout,st,en,thr,=sys.argv[1:]
table=pd.read_csv(filein,sep='\t')
st=int(float(st))
en=int(float(en))
thr=float(thr)
stcol=table.columns.to_list()[st]
encol=table.columns.to_list()[en]
table['length']=table[encol]-table[stcol]
subtable=table[table['length']>=thr]
subtable.to_csv(fileout,sep='\t',index=False,header=True)
