# importing
import pandas as pd
import sys
import itertools
import seaborn as sns
import matplotlib.pyplot as plt 

filein,fileout,lastev,randsize,ev1left,ev1right,ev2down,ev2up=sys.argv[1:]
lastev=int(float(lastev))
randsize=int(float(randsize))
ev1left=float(ev1left)
ev1right=float(ev1right)
ev2down=float(ev2down)
ev2up=float(ev2up)

# setting up
table=pd.read_csv(filein,sep='\t')
subtable=table[(table['EV1']>=ev1left)&(table['EV1']<=ev1right)]
subtable=subtable[(subtable['EV2']>=ev2down)&(subtable['EV2']<=ev2up)]
subtable=subtable.sample(n=randsize)
evs=list(table.columns[3:3+lastev])
pairs=itertools.combinations(evs,2)
tuplepairs=[pair for pair in pairs]
palette='tab20'

# plotting
clrs=['RACE','STUDY']
for clr in clrs:
	print(clr)
	for tup in tuplepairs:
		print(tup)
		fileoutseriously=fileout+'.'+clr+'.'+tup[0]+'.'+tup[1]+'.png'
		sns.scatterplot(subtable,x=tup[0],y=tup[1],hue=clr,palette=palette)
		plt.legend(title=clr)
		plt.savefig(fileoutseriously,dpi=300)
		plt.close()

