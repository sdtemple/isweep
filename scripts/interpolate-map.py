import pandas as pd
import argparse
import numpy as np
import scipy.interpolate as spi

# Set up argument parser
parser = argparse.ArgumentParser(description='Interpolate genetic map by cM step size.')
parser.add_argument('--map_input', 
                    type=str, 
                    help='Path to the input genetic map')
parser.add_argument('--file_output', 
                    type=str, 
                    help='Path to the output genetic map')
parser.add_argument('--cm_step_size', 
                    type=float,
                    default=0.02, 
                    help='cM increment for scan')

# Parse the arguments
args = parser.parse_args()
map_in = args.map_input
cm_step_size = args.cm_step_size
file_out = args.file_out
increment = float(cm_step_size)

mapping = pd.read_csv(map_in, header=None, sep='\t')
mapping.columns = ['chrom', 'rsid', 'cm', 'bp']

chromosomes = mapping['chrom'].unique()

submapping = mapping[mapping['chrom'] == chromosomes[0]]
smallest = submapping['cm'].min()
largest = submapping['cm'].max()

x = submapping['cm'].to_list()
y = submapping['bp'].to_list()

interp_func = spi.interp1d(x, y, kind='linear')
forecasts = np.arange(smallest + increment, largest - increment, increment)
weather = interp_func(forecasts)

new_dict = {
    'chrom': chromosomes[0],
    'rsid': '.',
    'cm': forecasts,
    'bp': weather
}
build_table = pd.DataFrame(new_dict)

for chromosome in chromosomes[1:]:

    # should not enter the for loop if only 1 chromosome in file

    submapping = mapping[mapping['chrom'] == chromosome]
    smallest = submapping['cm'].min()
    largest = submapping['cm'].max()

    x = submapping['cm'].to_list()
    y = submapping['bp'].to_list()

    interp_func = spi.interp1d(x, y, kind='linear')
    forecasts = np.arange(smallest + increment, largest - increment, increment)
    weather = interp_func(forecasts)

    new_dict = {
        'chrom': chromosome,
        'rsid': '.',
        'cm': forecasts,
        'bp': weather
    }
    new_table = pd.DataFrame(new_dict)

    build_table = pd.concat((build_table, new_table))

build_table.to_csv(file_out, sep='\t', index=False, header=False)