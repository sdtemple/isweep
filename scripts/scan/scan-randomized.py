# This file combines code from scripts/utilities/count-ibd-case.py and scripts/scan/multiple-testing-analytical-case.py
# To conduct many randomizations of the case/control labels and record the maximum value of the scan statistic
# And also to keep track of regions where the signal is often confounded (LCT in Europe, HLA in many populations)

# importing
import argparse
import numpy as np
import pandas as pd
from copy import deepcopy
from math import log, exp
import statsmodels.api as sm
import scipy
from scipy.stats import norm
from scipy.optimize import root_scalar
from scipy.stats import pearsonr
from random import shuffle

# Set up the argument parser
parser = argparse.ArgumentParser(description='Perform many scans with randomized case/control statuses.')
parser.add_argument('--output_scan_file', 
                    type=str,
                    required=True, 
                    help='Output scan data file path')
parser.add_argument('--output_sims_file', 
                    type=str,
                    required=True, 
                    help='Output max simulations file path')
parser.add_argument('--num_randomizations',
                    type=int,
                    required=True,
                    help='Number of randomizations to perform')
parser.add_argument('--input_case_file', 
                    type=str, 
                    required=True, 
                    help='File with case status for individuals')
parser.add_argument('--input_map_prefix', 
                    type=str,
                    required=True, 
                    help='Map file path prefix')
parser.add_argument('--input_map_suffix', 
                    type=str,
                    required=True, 
                    help='Map file path suffix')
parser.add_argument('--input_ibd_prefix', 
                    type=str,
                    required=True, 
                    help='Prefix of chromosome files')
parser.add_argument('--input_ibd_suffix', 
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
parser.add_argument('--Morgan_step_size', 
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
parser.add_argument('--start', 
                    type=int,
                    default=5, 
                    help='Start bp column index')
parser.add_argument('--end', 
                    type=int,
                    default=6, 
                    help='End bp column index')
parser.add_argument('--ind1', 
                    type=int, 
                    default=0, 
                    help='Individual 1 column index')
parser.add_argument('--ind2', 
                    type=int, 
                    default=2, 
                    help='Individual 2 column index')
parser.add_argument('--chunksize',
                    type=int,
                    default=10000000,
                    help='The parameter in pandas read_csv'
                    )

# Parse the arguments
args = parser.parse_args()
output_scan_file = args.output_scan_file
output_sims_file = args.output_sims_file
map_prefix = args.input_map_prefix
map_suffix = args.input_map_suffix
start = args.start
end = args.end
chunksize = args.chunksize
ind1 = args.ind1
ind2 = args.ind2
casefile = args.input_case_file
ibd_prefix = args.input_ibd_prefix
ibd_suffix = args.input_ibd_suffix
chr_st = args.chr_low
chr_en = args.chr_high
chr_len = args.chr_average_size
stepsize = args.Morgan_step_size
covariance_length = args.autocovariance_steps
pval = args.confidence_level
init_cut = args.outlier_cutoff
num_randomizations = args.num_randomizations

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

##### SETUP #####

# Initialize the scan data that we will keep
scan_data = {chr:None for chr in range(chr_st, chr_en + 1)}
copy_data = {chr:None for chr in range(chr_st, chr_en + 1)}
for chromosome in range(chr_st, chr_en + 1):
    try:

        # read in the map file
        map_file = pd.read_csv(f'{map_prefix}{chromosome}{map_suffix}', sep='\t', header=None)
        map_file.columns = ['chrom', 'rsid', 'cm', 'bp']

        # set up the scan data
        counts = np.zeros(map_file.shape[0],dtype=int)
        sigs = {
            'CHROM': np.ones(map_file.shape[0],dtype=int)*chromosome,
            'BPWINDOW': map_file['bp'].to_list(),
            'CMWINDOW': map_file['cm'].to_list(),
            'NUMSIG': counts,
        }
        scan_data[chromosome] = sigs

        # set up the copy data
        counts = np.zeros(map_file.shape[0],dtype=int)
        counter = {
            'CHROM': np.ones(map_file.shape[0],dtype=int)*chromosome,
            'BPWINDOW': map_file['bp'].to_list(),
            'CMWINDOW': map_file['cm'].to_list(),
            'COUNT0': deepcopy(counts),
            'COUNT1': deepcopy(counts),
        }
        copy_data[chromosome] = counter

    except:
        print('checkpoint 0')
        pass

#####

##### CASES AND CONTROLS #####

# Process the case file
case_df = pd.read_csv(casefile, sep='\t', header=None)
case_prop = case_df[1].mean()
sample_ids = case_df[0].to_list()
sample_num = len(sample_ids)
case_num = int(case_prop * sample_num)
pheno1 = [1 for _ in range(case_num)]
pheno0 = [0 for _ in range(sample_num - case_num)]
pheno = pheno1 + pheno0

#####

##### RANDOMIZATIONS #####

thrs_rand = []
thrs_anal = []
thrs_theta = []
for _ in range(num_randomizations):

    print(f'Randomization #{_+1} out of {num_randomizations}')

    # randomizing the case/control statuses
    shuffle(sample_ids)
    casedict = dict()
    for line in range(sample_num):
        casedict[str(sample_ids[line])] = pheno[line]

    copied_data = deepcopy(copy_data)
    final_goodbyes = []
    final_covs = []
    keep_all_data_control = []
    keep_all_data_case = []
    keep_bp = []
    keep_chr = []
    ctr = 0

    ##### COUNTING IBD SEGMENTS #####

    # counting case-case and control-control IBD segments
    for chromosome in range(chr_st, chr_en + 1):
        try:

            print('Starting chromosome ' + str(chromosome))

            # recall the basepair positions
            map_column = copied_data[chromosome]['BPWINDOW']

            # Counting IBD segments in chunks
            for chunk in pd.read_csv(f'{ibd_prefix}{chromosome}{ibd_suffix}', sep='\t', chunksize=chunksize, header=None):

                # Input and formatting
                columns = list(chunk.columns)
                encol = columns[end]
                stcol = columns[start]
                ind1col = columns[ind1]
                ind2col = columns[ind2]
                chunk[ind1col] = chunk[ind1col].astype(str)
                chunk[ind2col] = chunk[ind2col].astype(str)

                # Map case status to individuals
                chunk['case1'] = chunk[ind1col].map(casedict)
                chunk['case2'] = chunk[ind2col].map(casedict)
                chunk['match'] = chunk['case1'] == chunk['case2']
                chunk['casemult'] = chunk['case1'] * chunk['case2']
                chunk['case'] = chunk.apply(lambda x: x['casemult'] if x['match'] else 2, axis=1)

                # Counting IBD segments
                counts0 = [((chunk[stcol] <= i) & (chunk[encol] >= i) & (chunk['case'] == 0)).sum() for i in map_column]
                counts1 = [((chunk[stcol] <= i) & (chunk[encol] >= i) & (chunk['case'] == 1)).sum() for i in map_column]
                counts = np.array(counts)
                counts0 = np.array(counts0)
                counts1 = np.array(counts1)
                copied_data[chromosome]['COUNT0'] += counts0
                copied_data[chromosome]['COUNT1'] += counts1

            keep_all_data_control += copied_data[chromosome]['COUNT0'].tolist()
            keep_all_data_case += copied_data[chromosome]['COUNT1'].tolist()

        except Exception as e:
            print('Exception in chromosome ' + str(chromosome))
            pass

    #### PREPARE FOR NORMALIZATION #####

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

    ###### ESTIMATE THE COVARIANCE #####

    # Normalizing the data
    # And estimating the covariance
    # For each chromosome
    mx = -1e10
    for chromosome in range(chr_st, chr_en + 1):
        try:
            table = pd.DataFrame(copied_data[chromosome])
            subtable = table
            K = subtable.shape[0]
            ys1 = np.array(subtable['COUNT1'])
            ys0 = np.array(subtable['COUNT0'])
            zs1 = (ys1 - avg12) / std12 # normalize the two ibd rate processes
            zs0 = (ys0 - avg02) / std02
            xs = (zs1 - zs0 - avgk) / stdk # normalize the difference
            mx = max(mx, np.max(xs))

            byebye = covariance_length
            covs = []
            byes = []
            
            # loop over autoregressive lags
            for bye in range(1, byebye + 1):
                seq1 = np.arange(1, K - bye + 1, 1)
                seq2 = np.arange(1 + bye, K + 1, 1)
                # difference process
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

            final_goodbyes += byes
            final_covs += covs

            ctr += 1

        except:
            # this should not happen
            # this can happen when chromosome is too small
            # i.e. filtered out of the analysis
            print('Exception in chromosome ' + str(chromosome))
            pass
        
    final_goodbyes = np.array(final_goodbyes)
    final_covs = np.array(final_covs)
    # ignore the negatives b/c we will take the log
    # there should not be many estimated negative covariances
    # if there are, then the data is not well-behaved
    final_goodbyes = final_goodbyes[final_covs > 0.]
    final_covs = final_covs[final_covs > 0.]

    # difference process
    # fit the model [ minus log(cov) = theta * lag ]   
    xvar = np.array(final_goodbyes) * stepsize
    ytilde = - np.log(final_covs)
    fittilde = sm.OLS(ytilde, xvar).fit()
    theta = fittilde.params[0]

    # Calculating the multiple testing correction
    siglevel = pval
    chrnum = ctr
    chr_len = float(chr_len)
    lin2_root = root_scalar(OU_approx, args=(theta, stepsize, chrnum * chr_len, chrnum, siglevel), bracket=[1, 20])
    lin2 = lin2_root.root if lin2_root.converged else None
    if lin2 is not None:
        result_lin2 = lin2 * stdk + avgk
    else:
        # this should not happen
        print("Root not found for OU_approx")

    ##### KEEP TRACK OF THRESHOLDS #####

    thrs_anal.append(lin2)
    thrs_rand.append(mx)
    thrs_theta.append(theta)

    ##### KEEP TRACK OF REGIONS #####

    for chromosome in range(chr_st, chr_en + 1):
        try:
            table = pd.DataFrame(copied_data[chromosome])
            subtable = table
            K = subtable.shape[0]
            ys1 = np.array(subtable['COUNT1'])
            ys0 = np.array(subtable['COUNT0'])
            zs1 = (ys1 - avg12) / std12 # normalize the two ibd rate processes
            zs0 = (ys0 - avg02) / std02
            xs = (zs1 - zs0 - avgk) / stdk # normalize the difference
            bs = np.array(xs) >= lin2
            scan_data[chromosome]['NUMSIG'] += bs.astype(int)
        except:
            pass

    print('\n')

##### OUTPUT #####

with open(output_sims_file, 'w') as f:
    f.write('anal_z\trand_z\ttheta\n')
    for a, r, t in zip(thrs_anal, thrs_rand, thrs_theta):
        f.write(f'{a}\t{r}\t{t}\n')

output_table = pd.DataFrame(scan_data[chr_st])
for chromosome in range(chr_st + 1, chr_en + 1):
    try:
        output_table = pd.concat([output_table, pd.DataFrame(scan_data[chromosome])])
    except:
        pass
output_table.to_csv(output_scan_file, sep='\t', index=False)
