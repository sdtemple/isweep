import sys
import numpy as np

# these functions are provided by Bing AI (May 18, 2024)

def ecdf(data):
    x, counts = np.unique(data, return_counts=True)
    cumulative_sum = np.cumsum(counts)
    return x, cumulative_sum / cumulative_sum[-1]

def find_percentile(ecdf_x, ecdf_values, input_value):
    # Linear interpolation
    idx = np.searchsorted(ecdf_x, input_value)
    if idx == 0:
        return 0.0
    elif idx == len(ecdf_x):
        return 1.0
    else:
        x0, x1 = ecdf_x[idx - 1], ecdf_x[idx]
        y0, y1 = ecdf_values[idx - 1], ecdf_values[idx]
        slope = (y1 - y0) / (x1 - x0)
        return y0 + slope * (input_value - x0)

filein, qvalue, pvalue =sys.argv[1:]

qvalue = float(qvalue)
pvalue = float(pvalue)

rows = []
with open(filein,'r') as f:
    for line in f:
        row = [int(x) for x in line.strip().split('\t')]
        rows += row
rows = np.array(rows)

x,y = ecdf(rows)

out = find_percentile(x,y,qvalue)

cut1 = np.quantile(rows,out)
cut2 = np.quantile(rows,pvalue)
past1 = len(rows[rows>=cut1])
past2 = len(rows[rows>=cut2])
m = len(rows)

print(m)
print(out)
print(qvalue)
print(cut2)
print(past1)
print(past2)
print(past1/m)
print(past2/m)
print(1-pvalue)