# Multiple testing correction for Ornstein-Uhlenbeck approximation
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-09-24
# Description: This script is used to estimate a multiple testing correction,
# in an IBD case-control mapping, using the Ornstein-Uhlenbeck approximation.

# Import necessary modules
import numpy as np
import pandas as pd
from math import log, exp
import statsmodels.api as sm
import scipy
from scipy.stats import norm
from scipy.optimize import root_scalar
import argparse
from scipy.stats import pearsonr

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
parser.add_argument('file_out_test', 
                    type=str, 
                    help='Output file with multiple testing results')
parser.add_argument('file_out_cov', 
                    type=str, 
                    help='Output file with covariance estimates')
parser.add_argument('file_out_zdiff',
                    type=str,
                    help='Output file with test statistics')
parser.add_argument('file_out_cross_corr',
                    type=str,
                    help='Output file with cross-correlations per chromosome')
parser.add_argument('chr_prefix', 
                    type=str, 
                    help='Prefix of chromosome files')
parser.add_argument('chr_suffix', 
                    type=str, 
                    help='Suffix of chromosome files')
parser.add_argument('--chr_start', 
                    type=int, 
                    default=1, 
                    help='(default: 1) First chromosome number')
parser.add_argument('--chr_end', 
                    type=int, 
                    default=22, 
                    help='(default: 22) Last chromosome number')
parser.add_argument('--chr_len', 
                    type=float, 
                    default=1.5, 
                    help='(default: 1.5) Average length of chromosome (in Morgans)')
parser.add_argument('--chr_slide', 
                    type=float, 
                    default=0.0005, 
                    help='(default: 0.0005) Step size for each test (in Morgans)')
parser.add_argument('--cov_len', 
                    type=int, 
                    default=80, 
                    help='(default: 80) Number of steps for the covariance estimate')
parser.add_argument('--pvalue', 
                    type=float, 
                    default=0.05, 
                    help='(default: 0.05) p-value for the multiple testing correction')
parser.add_argument('--init_cut', 
                    type=float, 
                    default=4, 
                    help='(default: 4) Scalar for the initial outlier removal')
parser.add_argument('--counts_column0', 
                    type=str, 
                    default='COUNT0', 
                    help='(default: COUNT0) Column name of control-control counts')
parser.add_argument('--counts_column1', 
                    type=str, 
                    default='COUNT1', 
                    help='(default: COUNT1) Column name of case-case counts')
parser.add_argument('--distance_column',
                    type=str,
                    default='CMWINDOW',
                    help='(default: CMWINDOW) Column name of the distance between SNPs')

args = parser.parse_args()

# Assign arguments to variables
fileout = args.file_out_test
fileout2 = args.file_out_cov
prefix = args.chr_prefix
suffix = args.chr_suffix
st = args.chr_start
en = args.chr_end
chr_len = args.chr_len
stepsize = args.chr_slide
covariance_length = args.cov_len
pval = args.pvalue
init_cut = args.init_cut
counts_column0 = args.counts_column0
counts_column1 = args.counts_column1

final_goodbyes = []
final_covs = []
keep_all_data_control = []
keep_all_data_case = []
keep_distance = []
keep_chr = []
ctr = 0

# First pass to determine the mean and standard deviation
for chromosome in range(st, en + 1):
    try:
        table = pd.read_csv(prefix + str(chromosome) + suffix, sep='\t')
        subtable = table
        keep_all_data_control += subtable[counts_column0].to_list()
        keep_all_data_case += subtable[counts_column1].to_list()
        keep_chr += [chromosome] * subtable.shape[0]
        keep_distance += subtable[args.distance_column].to_list()
    except:
        pass
keep_all_data_case = np.array(keep_all_data_case)
keep_all_data_control = np.array(keep_all_data_control)
avg1 = np.median(keep_all_data_case) # median more robust, but mean can be used
std1 = np.std(keep_all_data_case)
avg0 = np.median(keep_all_data_control) # expect that there are outliers (putatively selection)
std0 = np.std(keep_all_data_control)

# Remove obvious outliers
upper0 = avg0 + init_cut * std0
lower0 = avg0 - init_cut * std0
upper1 = avg1 + init_cut * std1
lower1 = avg1 - init_cut * std1
keep_all_data_case2 = keep_all_data_case[(abs(keep_all_data_case - avg1) < (init_cut * std1)) * (keep_all_data_case > 0)]
keep_all_data_control2 = keep_all_data_control[(abs(keep_all_data_control - avg0) < (init_cut * std0)) * (keep_all_data_control > 0)]

# Recalculate the mean and standard deviation
std12 = np.std(keep_all_data_case2)
avg12 = np.mean(keep_all_data_case2) # use the mean now
std02 = np.std(keep_all_data_control2)
avg02 = np.mean(keep_all_data_control2)
# will use avg12, std12, avg02, std02 to normalize later

# normalize the two ibd rate processes
k1 = (keep_all_data_case - avg12) / std12
k0 = (keep_all_data_control - avg02) / std02

# take the difference
kz = k1 - k0

# prepare to output a normalized version of the test statistic
out = dict()
out['CHROM'] = keep_chr
out['DISTANCE'] = keep_distance
out['Z0'] = k0
out['Z1'] = k1 
out['ZDIFF'] = kz

# compute mean and standard deviation of the difference
kz = kz[(abs(keep_all_data_case - avg1) < init_cut * std1)
        * (keep_all_data_case > 0)
        * (abs(keep_all_data_control - avg0) < init_cut * std0)
        * (keep_all_data_control > 0)
        ]
avgk = np.mean(kz) # mean should be zero by design
stdk = np.std(kz)
out['ZDIFFZ'] = (k1 - k0 - avgk) / stdk
# will use avgk and stdk to normalize later

# output normalized test statistic
outpd = pd.DataFrame(out)
outpd.to_csv(args.file_out_zdiff, sep='\t', index=False)

# Open an output file
# And write the header
# Which is the time lags
g = open(fileout2, 'w')
for c in range(1, covariance_length + 1):
    g.write(str(c * stepsize))
    g.write('\t')
g.write('\n')

cc = open(args.file_out_cross_corr, 'w')
cc.write('CHROM\tCROSS_CORR\n')

# Normalizing the data
# And estimating the covariance
# For each chromosome
cross = []
for chromosome in range(st, en + 1):
    try:
        table = pd.read_csv(prefix + str(chromosome) + suffix, sep='\t')
        subtable = table
        K = subtable.shape[0]
        ys1 = np.array(subtable[counts_column1])
        ys0 = np.array(subtable[counts_column0])
        zs1 = (ys1 - avg12) / std12 # normalize the two ibd rate processes
        zs0 = (ys0 - avg02) / std02
        xs = (zs1 - zs0 - avgk) / stdk # normalize the difference

        # calculate the cross correlation
        zs1c = zs1[(abs(ys1 - avg1) < init_cut * std1) * (ys1 > 0) * (abs(ys0 - avg0) < init_cut * std0) * (ys0 > 0)]
        zs0c = zs0[(abs(ys1 - avg1) < init_cut * std1) * (ys1 > 0) * (abs(ys0 - avg0) < init_cut * std0) * (ys0 > 0)]
        cross = pearsonr(zs1c, zs0c)[0]
        cc.write(str(chromosome) + '\t' + str(cross) + '\n')

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
            bools = [ys1[x]<upper1 
                     and ys1[x]>lower1 
                     and ys1[x]>0 
                     and ys0[x]<upper0 
                     and ys0[x]>lower0 
                     and ys0[x]>0
                     and ys1[x+bye]<upper1
                     and ys1[x+bye]>lower1
                     and ys1[x+bye]>0
                     and ys0[x+bye]<upper0
                     and ys0[x+bye]>lower0
                     and ys0[x+bye]>0 
                     for x in range(ln)
                     ]
            xs1 = xs1[bools]
            xs2 = xs2[bools]
            # calculate the covariance
            covs.append(np.cov(xs1, xs2)[0, 1])
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
    result_lin = lin * stdk + avgk
else:
    # this should not happen
    print("Root not found for OU_approx_cont")

if lin2 is not None:
    result_lin2 = lin2 * stdk + avgk
else:
    # this should not happen
    print("Root not found for OU_approx")

# Write the results to an output file

f = open(fileout,'w')
f.write('p-value:'); f.write(str(pval)); f.write('\n')
f.write('chromosome-number:'); f.write(str(chrnum)); f.write('\n')
f.write('average-chromosome-length-morgan:'); f.write(str(chr_len)); f.write('\n')
f.write('step-size-morgan:'); f.write(str(stepsize)); f.write('\n')
f.write('estimated-theta:'); f.write(str(theta)); f.write('\n')
f.write('upper-discrete-z:'); f.write(str(lin2)); f.write('\n')
f.write('upper-continuous-z:'); f.write(str(lin)); f.write('\n')
f.write('upper-discrete-raw:'); f.write(str(result_lin2)); f.write('\n')
f.write('upper-continuous-raw:'); f.write(str(result_lin)); f.write('\n')
f.write('case-lower-bound:'); f.write(str(avg1 - init_cut * std1)); f.write('\n')
f.write('case-upper-bound:'); f.write(str(avg1 + init_cut * std1)); f.write('\n')
f.write('control-lower-bound:'); f.write(str(avg0 - init_cut * std0)); f.write('\n')
f.write('control-upper-bound:'); f.write(str(avg0 + init_cut * std0)); f.write('\n')
f.write('case-revised-mean:'); f.write(str(avg12)); f.write('\n')
f.write('case-revised-std:'); f.write(str(std12)); f.write('\n')
f.write('control-revised-mean:'); f.write(str(avg02)); f.write('\n')
f.write('control-revised-std:'); f.write(str(std02)); f.write('\n')
kza = np.array(k1-k0)
# kza = abs(kza)
kzz = np.array((k1 - k0 - avgk) / stdk)
# kzz = abs(kzz)
f.write('actual-maximum-z:'); f.write(str(np.max(kzz))); f.write('\n')
f.write('actual-maximum-raw:'); f.write(str(np.max(kza))); f.write('\n')
f.close()
g.close()
cc.close()
