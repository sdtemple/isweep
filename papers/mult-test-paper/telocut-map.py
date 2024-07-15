import pandas as pd
import sys
map_in,map_out,telocut=sys.argv[1:]
table=pd.read_csv(map_in,sep='\t')
cmcol=table.columns.to_list()[2]
cut=float(telocut)
mn=table[cmcol].min()+cut
mx=table[cmcol].max()-cut
subtable=table[(table[cmcol]>=mn)&(table[cmcol]<=mx)]
subtable.to_csv(map_out,sep='\t',index=False,header=True)