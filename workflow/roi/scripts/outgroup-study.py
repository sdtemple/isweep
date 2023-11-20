# setup
import sys
import numpy as np
import pandas as pd
folder, filein = sys.argv[1:]

# read outlier files
ctr= 0 
stu = pd.read_csv(filein, sep='\t')

try:
    while True:
        ctr += 1
        print('outlier', '', ctr)
        tab = pd.read_csv(folder + '/outlier' + str(ctr) + '.txt', sep='\t')
        tab.columns = ['SAMPLE_ID']    
        tab['SAMPLE_ID'] = tab['SAMPLE_ID'].apply(lambda x: str(x)[:-2])
        tab = tab.merge(stu, on='SAMPLE_ID')
        print(tab['STUDY'].value_counts())
        print('\n')
except:
    pass