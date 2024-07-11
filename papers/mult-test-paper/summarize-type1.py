import sys
import numpy as np
import pandas as pd
fileout,prefix,suffix,k,j,=sys.argv[1:]
k=int(float(k))
j=int(float(j))
dat=dict()
for i in range(k,j+1):
    try:
        with open(prefix+str(i)+suffix,'r') as f:
            for line in f:
                Key,Value=line.strip().split(':')
                # print(Key)
                try:
                    dat[Key].append(float(Value))
                except:
                    dat[Key] = [float(Value)]
    except:
        print(i)
actual=np.array(dat['actual-max'])
disc=np.array(dat['upper-discrete-raw'])
cont=np.array(dat['upper-continuous-raw'])
table=pd.DataFrame(dat)
table.to_csv(fileout,sep='\t',index=False,header=True)
print(np.mean(actual>=disc))
