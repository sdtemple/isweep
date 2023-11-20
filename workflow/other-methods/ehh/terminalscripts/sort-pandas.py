import pandas as pd
import sys
tablein,tableout=sys.argv[1:]
tab=pd.read_csv(tablein,sep='\t')
tab.sort_values(by='iSAFE',ascending=False,inplace=True)
tab.reset_index(inplace=True)
tab['ZDELTA']=tab['iSAFE']
tab['AAF']=tab['DAF']
tab.to_csv(tableout,sep='\t',index=False)
