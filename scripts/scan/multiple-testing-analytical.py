# Multiple testing correction for Ornstein-Uhlenbeck approximation
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-07-18
# Description: This script is used to estimate a multiple testing correction,
# in an IBD-based selection scan, using the Ornstein-Uhlenbeck approximation.

# Import necessary modules
import numpy as np
import pandas as pd
from math import log, exp
import statsmodels.api as sm
import scipy
from scipy.stats import norm
from scipy.optimize import root_scalar
import argparse

# Functions for the OU approximation (Ornstein-Uhlenbeck)
def Nu(y):
    y = y / 2
    return (1 / y) * (norm.cdf(y) - 0.5) / (y * norm.cdf(y) + norm.pdf(y))

def OU_approx(z, beta, Delta, length, chr, center, test="one-sided"):
    d = 2 if test == "two-sided" else 1
    p = 1 - np.exp(-d * chr * (1 - norm.cdf(z)) - d * beta * length * z * norm.pdf(z) * Nu(z * np.sqrt(2 * beta * Delta)))
    return p - center

def OU_approx_cont(z, beta, length, chr, center, test="one-sided"):
    d = 2 if test == "two-sided" else 1
    p = 1 - np.exp(-d * chr * (1 - norm.cdf(z)) - d * beta * length * z * norm.pdf(z))
    return p - center

# Set up argument parser
parser = argparse.ArgumentParser(description="Estimate a multiple testing correction using the Ornstein-Uhlenbeck approximation.")
parser.add_argument('--output_testing_file', 
                    type=str,
                    required=True, 
                    help='Output file with multiple testing results')
parser.add_argument('--output_autocov_file', 
                    type=str,
                    required=True,
                    help='Output file with covariance estimates')
parser.add_argument('--input_prefix', 
                    type=str,
                    required=True, 
                    help='Prefix of chromosome files')
parser.add_argument('--input_suffix', 
                    type=str,
                    required=True, 
                    help='Suffix of chromosome files')
parser.add_argument('--chr_low', 
                    type=int, 
                    default=1, 
                    help='(default: 1) First chromosome number')
parser.add_argument('--chr_high', 
                    type=int, 
                    default=22, 
                    help='(default: 22) Last chromosome number')
parser.add_argument('--chr_average_size', 
                    type=float, 
                    default=1.5, 
                    help='(default: 1.5) Average length of chromosome (in Morgans)')
parser.add_argument('--cM_step_size', 
                    type=float, 
                    default=0.0005, 
                    help='(default: 0.0005) Step size for each test (in Morgans)')
parser.add_argument('--autocovariance_steps', 
                    type=int, 
                    default=80, 
                    help='(default: 80) Number of steps for the covariance estimate')
parser.add_argument('--confidence_level', 
                    type=float, 
                    default=0.05, 
                    help='(default: 0.05) p-value for the multiple testing correction')
parser.add_argument('--outlier_cutoff', 
                    type=float, 
                    default=4, 
                    help='(default: 4) Scalar for the initial outlier removal')
parser.add_argument('--counts_column', 
                    type=str, 
                    default='COUNT', 
                    help='(default: COUNT) Column name of counts')

args = parser.parse_args()

# Assign arguments to variables
fileout = args.output_testing_file
fileout2 = args.output_autocov_file
prefix = args.input_prefix
suffix = args.input_suffix
st = args.chr_low
en = args.chr_high
chr_len = args.chr_average_size
stepsize = args.cM_step_size
covariance_length = args.autocovariance_steps
pval = args.confidence_level
init_cut = args.outlier_cutoff
counts_column = args.counts_column

final_goodbyes = []
final_covs = []
keep_all_data = []
ctr = 0

# First pass to determine the mean and standard deviation
for chromosome in range(st, en + 1):
    try:
        table = pd.read_csv(prefix + str(chromosome) + suffix, sep='\t')
        subtable = table
        keep_all_data += subtable[counts_column].to_list()
    except:
        pass
stddev = np.std(keep_all_data)
# avg = np.mean(keep_all_data)
avg = np.median(keep_all_data)
print('first pass: sigma, mu, upper, lower')
print(stddev)
print(avg)
print(avg + init_cut * stddev)
print(avg - init_cut * stddev)

# Second pass to determine the mean and standard deviation
# After removing obvious outliers
keep_all_data2 = [k for k in keep_all_data if k < avg + init_cut * stddev and k > avg - init_cut * stddev and k > 0]
stddev2 = np.std(keep_all_data2)
avg2 = np.mean(keep_all_data2)
print('\n\n\n')
print('second pass: sigma, mu, upper, lower')
print(stddev2)
print(avg2)
print(avg2 + init_cut * stddev2)
print(avg2 - init_cut * stddev2)
print('\n\n\n')

# This is the upper bound for the initial outliers
upper = avg + init_cut * stddev
lower = avg - init_cut * stddev

# Open an output file
# And write the header
# Which is the time lags
g = open(fileout2, 'w')
for c in range(1, covariance_length + 1):
    g.write(str(c * stepsize))
    g.write('\t')
g.write('\n')

# Normalizing the data
# And estimating the covariance
# For each chromosome
for chromosome in range(st, en + 1):
    try:
        table = pd.read_csv(prefix + str(chromosome) + suffix, sep='\t')
        subtable = table
        K = subtable.shape[0]
        xs = np.array(subtable[counts_column])

        byebye = covariance_length
        covs = []
        byes = []
        
        # loop over autoregressive lags
        for bye in range(1, byebye + 1):
            seq1 = np.arange(1, K - bye + 1, 1)
            seq2 = np.arange(1 + bye, K + 1, 1)
            xs1 = xs[seq1-1]  # Adjusting index to 0-based for Python
            xs2 = xs[seq2-1]
            ln = len(xs1)
            # ignore the outliers
            bools = [xs1[x] < upper and xs1[x] > lower and xs1[x] > 0 and xs2[x] < upper and xs2[x] > lower and xs2[x] > 0 for x in range(ln)]
            xs1 = xs1[bools]
            xs2 = xs2[bools]
            # normalizing
            zs1 = (xs1 - avg2) / stddev2
            zs2 = (xs2 - avg2) / stddev2
            # calculate the covariance
            covs.append(np.cov(zs1, zs2)[0, 1])
            byes.append(bye)

        # write the covariance estimates to an output file
        for c in covs:
            g.write(str(c))
            g.write('\t')
        g.write('\n')

        final_goodbyes += byes
        final_covs += covs

        ctr += 1

    except:
        # this should not happen
        print('here')
        pass
    
final_goodbyes = np.array(final_goodbyes)
final_covs = np.array(final_covs)
# ignore the negatives b/c we will take the log
# there should not be many estimated negative covariances
# if there are, then the data is not well-behaved
final_goodbyes = final_goodbyes[final_covs > 0.]
final_covs = final_covs[final_covs > 0.]

# fit the model [ minus log(cov) = theta * lag ]   
xvar = np.array(final_goodbyes) * stepsize
ytilde = - np.log(final_covs)
fittilde = sm.OLS(ytilde, xvar).fit()

print(fittilde.summary())
print('\n\n\n\n')

print('theta')
print(fittilde.params[0])
print('\n\n\n\n')

# Calculating the multiple testing correction

siglevel = pval
chrnum = ctr
chr_len = float(chr_len)
theta = fittilde.params[0]

# Using scipy's root_scalar to find the roots similar to R's uniroot
lin_root = root_scalar(OU_approx_cont, args=(theta, chrnum * chr_len, chrnum, siglevel), bracket=[1, 20])
lin = lin_root.root if lin_root.converged else None

lin2_root = root_scalar(OU_approx, args=(theta, stepsize, chrnum * chr_len, chrnum, siglevel), bracket=[1, 20])
lin2 = lin2_root.root if lin2_root.converged else None

if lin is not None:
    result_lin = lin * stddev2 + avg2
else:
    # this should not happen
    print("Root not found for OU_approx_cont")

if lin2 is not None:
    result_lin2 = lin2 * stddev2 + avg2
else:
    # this should not happen
    print("Root not found for OU_approx")

# Write the results to an output file

f = open(fileout,'w')
f.write('confidence-level:\t'); f.write(str(pval)); f.write('\n')
f.write('chromosome-number:\t'); f.write(str(chrnum)); f.write('\n')
f.write('average-chromosome-length-morgan:\t'); f.write(str(chr_len)); f.write('\n')
f.write('step-size-morgan:\t'); f.write(str(stepsize)); f.write('\n')
f.write('estimated-theta:\t'); f.write(str(theta)); f.write('\n')
f.write('revised-mean:\t'); f.write(str(avg2)); f.write('\n')
f.write('revised-standard-deviation:\t'); f.write(str(stddev2)); f.write('\n')
f.write('upper-discrete-z:\t'); f.write(str(lin2)); f.write('\n')
f.write('upper-continuous-z:\t'); f.write(str(lin)); f.write('\n')
f.write('upper-discrete-raw:\t'); f.write(str(result_lin2)); f.write('\n')
f.write('upper-continuous-raw:\t'); f.write(str(result_lin)); f.write('\n')
f.write('initial-mean:\t'); f.write(str(avg)); f.write('\n')
f.write('initial-standard-deviation:\t'); f.write(str(stddev)); f.write('\n')
f.write('initial-scalar:\t'); f.write(str(init_cut)); f.write('\n')
f.write('initial-upper-bound:\t'); f.write(str(upper)); f.write('\n')
f.write('initial-lower-bound:\t'); f.write(str(lower)); f.write('\n')
f.close()
g.close()
