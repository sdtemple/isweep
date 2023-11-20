import gzip
import sys
filein,fileout=sys.argv[1:]
g = gzip.open(fileout, 'wt')
f = gzip.open(filein, 'rt')
g.write(f.readline())
for line in f:
    cells = line.strip().split('\t')
    if cells[-1] != '>3':
        g.write(line)
g.close()
f.close()

# ### probably slower

# import sys
# import pandas as pd

# filein,fileout=sys.argv[1:]
# table=pd.read_csv(filein,sep='\t')
# table['degree2']=table['degree'].replace('>3',4).astype('float')
# subtable=table[table['degree2']<4.]
# subtable.to_csv(fileout,sep='\t',index=False,header=True)



        
        