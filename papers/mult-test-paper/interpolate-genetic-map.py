import pandas as pd
import sys
import numpy as np
import scipy.interpolate as spi

maps, cm_slide, fileout = sys.argv[1:]
increment = float(cm_slide)

mapping = pd.read_csv(maps,header=None,sep='\t')
mapping.columns = ['chrom','rsid','cm','bp']
print(mapping.head())

chromosomes = mapping['chrom'].unique()
print(chromosomes)

submapping = mapping[mapping['chrom']==chromosomes[0]]
smallest = submapping['cm'].min()
largest = submapping['cm'].max()

x = submapping['cm'].to_list()
y = submapping['bp'].to_list()

interp_func = spi.interp1d(x, y, kind='linear')
forecasts = np.arange(smallest+increment,largest-increment,increment)
weather = interp_func(forecasts)

new_dict = dict()
new_dict['chrom'] = chromosomes[0]
new_dict['rsid'] = '.'
new_dict['cm'] = forecasts
new_dict['bp'] = weather
build_table = pd.DataFrame(new_dict)

for chromosome in chromosomes[1:]:

    submapping = mapping[mapping['chrom']==chromosome]
    smallest = submapping['cm'].min()
    largest = submapping['cm'].max()

    x = submapping['cm'].to_list()
    y = submapping['bp'].to_list()

    interp_func = spi.interp1d(x, y, kind='linear')
    forecasts = np.arange(smallest+increment,largest-increment,increment)
    weather = interp_func(forecasts)

    new_dict = dict()
    new_dict['chrom'] = chromosome
    new_dict['rsid'] = '.'
    new_dict['cm'] = forecasts
    new_dict['bp'] = weather
    new_table = pd.DataFrame(new_dict)

    build_table = pd.concat((build_table,new_table))

build_table.to_csv(fileout,sep='\t',index=False)