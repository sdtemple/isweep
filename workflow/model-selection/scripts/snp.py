import sys
import pandas as pd
tablein,fileout,low=sys.argv[1:]
low=float(low)
table=pd.read_csv(tablein,sep='\t')
table=table[table['AAF']>=low]
bestbp=table['POS'].tolist()[0]
bestaf=table['AAF'].tolist()[0]
f=open(fileout,'w')
f.write('bp\t')
f.write(str(int(bestbp)))
f.write('\n')
f.write('frequency\t')
f.write(str(round(bestaf,4)))
f.write('\n')
f.close()
