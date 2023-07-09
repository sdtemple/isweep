import sys
import pandas as pd
tablein,fileout=sys.argv[1:]
table=pd.read_csv(tablein,sep='\t')
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
