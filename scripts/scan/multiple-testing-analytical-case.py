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
parser = argparse.ArgumentParser(description="Estimate a multiple testing correction using the Ornstein-Uhlenbeck approximation (case-control mapping).")

# Set up argument parser
parser = argparse.ArgumentParser(description="Estimate a multiple testing correction using the Ornstein-Uhlenbeck approximation (case-control mapping).")
parser.add_argument('--output_testing_file', 
                    type=str,
                    required=True, 
                    help='Output file with multiple testing results')
parser.add_argument('--output_autocov_prefix', 
                    type=str,
                    required=True,
                    help='Output file prefix with autocovariance estimates')
parser.add_argument('--output_crosscov_file', 
                    type=str,
                    required=True, 
                    help='Output file with crosscovariance estimates')
parser.add_argument('--output_scan_file', 
                    type=str,
                    required=True, 
                    help='Output file standardized scanning statistics')
parser.add_argument('--output_excess_file', 
                    type=str,
                    required=True, 
                    help='Output file standardized scanning statistics if GW significant')
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
parser.add_argument('--counts_column0', 
                    type=str, 
                    default='COUNT0', 
                    help='(default: COUNT0) Column name of control-control counts')
parser.add_argument('--counts_column1', 
                    type=str, 
                    default='COUNT1', 
                    help='(default: COUNT1) Column name of case-case counts')

args = parser.parse_args()

# Assign arguments to variables
fileout = args.output_testing_file
fileout2 = args.output_autocov_prefix
fileout3 = fileout2 + '.case.tsv'
fileout6 = fileout2 + '.control.tsv'
fileout2 = fileout2 + '.diff.tsv'
prefix = args.input_prefix
suffix = args.input_suffix
st = args.chr_low
en = args.chr_high
chr_len = args.chr_average_size
stepsize = args.cM_step_size
covariance_length = args.autocovariance_steps
pval = args.confidence_level
init_cut = args.outlier_cutoff
counts_column0 = args.counts_column0
counts_column1 = args.counts_column1

final_goodbyes = []
final_covs = []
final_covs0 = []
final_covs1 = []
keep_all_data_control = []
keep_all_data_case = []
keep_distance = []
keep_bp = []
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
        keep_distance += subtable['CMWINDOW'].to_list()
        keep_bp += subtable['BPWINDOW'].to_list()
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

# compute mean and standard deviation of the difference
kz = kz[(abs(keep_all_data_case - avg1) < init_cut * std1)
        * (keep_all_data_case > 0)
        * (abs(keep_all_data_control - avg0) < init_cut * std0)
        * (keep_all_data_control > 0)
        ]
avgk = np.mean(kz) # mean should be zero by design
stdk = np.std(kz)
# will use avgk and stdk to normalize later

# take the difference
kz = k1 - k0

# prepare to output a normalized version of the test statistic
out = dict()
out['BPWINDOW'] = keep_bp
out['CMWINDOW'] = keep_distance
out['ZDIFFZ'] = (k1 - k0 - avgk) / stdk
out['CHROM'] = keep_chr
out['Z0'] = k0
out['Z1'] = k1 
out['COUNT0'] = keep_all_data_control
out['COUNT1'] = keep_all_data_case
out['ZDIFF'] = kz

# difference process
# Open an output file
# And write the header
# Which is the time lags
g = open(fileout2, 'w')
for c in range(1, covariance_length + 1):
    g.write(str(c * stepsize))
    g.write('\t')
g.write('\n')


# control process
# Open an output file
# And write the header
# Which is the time lags
g0 = open(fileout6, 'w')
for c in range(1, covariance_length + 1):
    g0.write(str(c * stepsize))
    g0.write('\t')
g0.write('\n')

# case process
# Open an output file
# And write the header
# Which is the time lags
g1 = open(fileout3, 'w')
for c in range(1, covariance_length + 1):
    g1.write(str(c * stepsize))
    g1.write('\t')
g1.write('\n')

cc = open(args.output_crosscov_file, 'w')
cc.write('CHROM\tCROSS_COV\n')

# Normalizing the data
# And estimating the covariance
# For each chromosome
cross = []
cross_case = []
cross_control = []
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
        cross_case += list(zs1c)
        cross_control += list(zs0c)
        cross = pearsonr(zs1c, zs0c)[0]
        cc.write(str(chromosome) + '\t' + str(cross) + '\n')

        byebye = covariance_length
        covs = []
        covs0 = []
        covs1 = []
        byes = []
        
        # loop over autoregressive lags
        for bye in range(1, byebye + 1):
            seq1 = np.arange(1, K - bye + 1, 1)
            seq2 = np.arange(1 + bye, K + 1, 1)
            # difference process
            xs1 = xs[seq1-1]  # Adjusting index to 0-based for Python
            xs2 = xs[seq2-1]
            # control process
            us1 = zs0[seq1-1]  # Adjusting index to 0-based for Python
            us2 = zs0[seq2-1]
            # case process
            vs1 = zs1[seq1-1]  # Adjusting index to 0-based for Python
            vs2 = zs1[seq2-1]
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
            us1 = us1[bools]
            us2 = us2[bools]
            vs1 = vs1[bools]
            vs2 = vs2[bools]
            # calculate the covariance
            covs.append(np.cov(xs1, xs2)[0, 1])
            covs0.append(np.cov(us1, us2)[0, 1])
            covs1.append(np.cov(vs1, vs2)[0, 1])
            byes.append(bye)

        # write the covariance estimates to an output file
        for c in covs:
            g.write(str(c))
            g.write('\t')
        g.write('\n')

        # write the covariance estimates to an output file
        for c in covs0:
            g1.write(str(c))
            g1.write('\t')
        g1.write('\n')

        # write the covariance estimates to an output file
        for c in covs1:
            g0.write(str(c))
            g0.write('\t')
        g0.write('\n')

        final_goodbyes += byes
        final_covs += covs

        final_covs0 += covs0
        final_covs1 += covs1

        ctr += 1

    except:
        # this should not happen
        print('here')
        pass
    
final_goodbyes = np.array(final_goodbyes)
final_covs = np.array(final_covs)
final_covs0 = np.array(final_covs0)
final_covs1 = np.array(final_covs1)
# ignore the negatives b/c we will take the log
# there should not be many estimated negative covariances
# if there are, then the data is not well-behaved
final_goodbyes0 = final_goodbyes[final_covs0 > 0.]
final_goodbyes1 = final_goodbyes[final_covs1 > 0.]
final_covs0 = final_covs0[final_covs0 > 0.]
final_covs1 = final_covs1[final_covs1 > 0.]
final_goodbyes = final_goodbyes[final_covs > 0.]
final_covs = final_covs[final_covs > 0.]

cross_corr_all, _ = pearsonr(np.array(cross_case), np.array(cross_control))

# difference process
# fit the model [ minus log(cov) = theta * lag ]   
xvar = np.array(final_goodbyes) * stepsize
ytilde = - np.log(final_covs)
fittilde = sm.OLS(ytilde, xvar).fit()

print(fittilde.summary())
print('\n\n\n\n')

print('theta')
print(fittilde.params[0])
print('\n\n\n\n')

theta = fittilde.params[0]

# control process
# fit the model [ minus log(cov) = theta * lag ]   
xvar = np.array(final_goodbyes0) * stepsize
ytilde = - np.log(final_covs0)
fittilde = sm.OLS(ytilde, xvar).fit()

print(fittilde.summary())
print('\n\n\n\n')

print('theta')
print(fittilde.params[0])
print('\n\n\n\n')

theta0 = fittilde.params[0]

# case process
# fit the model [ minus log(cov) = theta * lag ]   
xvar = np.array(final_goodbyes1) * stepsize
ytilde = - np.log(final_covs1)
fittilde = sm.OLS(ytilde, xvar).fit()

print(fittilde.summary())
print('\n\n\n\n')

print('theta')
print(fittilde.params[0])
print('\n\n\n\n')

theta1 = fittilde.params[0]

# Calculating the multiple testing correction

siglevel = pval
chrnum = ctr
chr_len = float(chr_len)

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
f.write('confidence-level:\t'); f.write(str(pval)); f.write('\n')
f.write('chromosome-number:\t'); f.write(str(chrnum)); f.write('\n')
f.write('average-chromosome-length-morgan:\t'); f.write(str(chr_len)); f.write('\n')
f.write('step-size-morgan:\t'); f.write(str(stepsize)); f.write('\n')
f.write('estimated-theta:\t'); f.write(str(theta)); f.write('\n')
f.write('estimated-theta0:\t'); f.write(str(theta0)); f.write('\n')
f.write('estimated-theta1:\t'); f.write(str(theta1)); f.write('\n')
f.write('estimated-rho:\t'); f.write(str(cross_corr_all)); f.write('\n')
f.write('case-revised-mean:\t'); f.write(str(avg12)); f.write('\n')
f.write('case-revised-std:\t'); f.write(str(std12)); f.write('\n')
f.write('control-revised-mean:\t'); f.write(str(avg02)); f.write('\n')
f.write('control-revised-std:\t'); f.write(str(std02)); f.write('\n')
f.write('upper-discrete-z:\t'); f.write(str(lin2)); f.write('\n')
f.write('upper-continuous-z:\t'); f.write(str(lin)); f.write('\n')
f.write('upper-discrete-raw:\t'); f.write(str(result_lin2)); f.write('\n')
f.write('upper-continuous-raw:\t'); f.write(str(result_lin)); f.write('\n')
f.write('initial-case-lower-bound:\t'); f.write(str(avg1 - init_cut * std1)); f.write('\n')
f.write('initial-case-upper-bound:\t'); f.write(str(avg1 + init_cut * std1)); f.write('\n')
f.write('initial-control-lower-bound:\t'); f.write(str(avg0 - init_cut * std0)); f.write('\n')
f.write('initial-control-upper-bound:\t'); f.write(str(avg0 + init_cut * std0)); f.write('\n')
kza = np.array(k1-k0)
# kza = abs(kza)
kzz = np.array((k1 - k0 - avgk) / stdk)
# kzz = abs(kzz)
f.write('actual-maximum-z:\t'); f.write(str(np.max(kzz))); f.write('\n')
f.write('actual-maximum-raw:\t'); f.write(str(np.max(kza))); f.write('\n')
f.close()
g.close()
g1.close()
g0.close()
cc.close()

# output normalized test statistic with thresholds
out['PVALUE'] = norm.sf(out['ZDIFFZ'])
out['UPPER_ANALYTICAL'] = lin2
out['Z_UPPER_ANALYTICAL'] = result_lin2
out['GW_LEVEL_ANALYTICAL'] = norm.sf(result_lin2)
out['ADJ_CASE_MEAN'] = avg12
out['ADJ_CASE_STDDEV'] = std12
out['ADJ_CONTROL_MEAN'] = avg02
out['ADJ_CONTROL_STDDEV'] = std02
out['INIT_CASE_LOWER_BOUND'] = avg1 - init_cut * std1
out['INIT_CASE_UPPER_BOUND'] = avg1 + init_cut * std1
out['INIT_CONTROL_LOWER_BOUND'] = avg0 - init_cut * std0
out['INIT_CONTROL_UPPER_BOUND'] = avg0 + init_cut * std0
out['UPPER_CONTINUOUS'] = lin
out['Z_UPPER_CONTINUOUS'] = result_lin
out['GW_LEVEL_CONTINUOUS'] = norm.sf(result_lin)
out['CONFLEVEL'] = pval

outpd = pd.DataFrame(out)
outpd.sort_values(['CHROM','CMWINDOW'],inplace=True)
outpd.to_csv(args.output_scan_file, sep='\t', index=False)

# output the significant regions
# Filtered subtable
suboutpd = outpd[outpd['ZDIFFZ'] > lin2]
suboutpd.to_csv(args.output_excess_file, sep='\t', header=True, index=False)
