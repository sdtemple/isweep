import pandas as pd
import warnings
import sys
warnings.filterwarnings('ignore')

# inputs
filein1,filein2,fileout=sys.argv[1:]
table1=pd.read_csv(filein1,header=None,sep='\t')
table2=pd.read_csv(filein2,header=None,sep='\t')
table1.columns=['CHROM','POS']
table2.columns=['CHROM','POS']

# outputs
assert table1['CHROM'][0] == table2['CHROM'][0]
assert len(table1['CHROM'].unique()) == 1
assert len(table2['CHROM'].unique()) == 1
subtable=table1[table1['POS'].isin(table2['POS'])]
subtable=subtable[['CHROM','POS']]
subtable.columns=['CHROM','POS']
subtable.to_csv(fileout,sep='\t',header=False,index=False)

