import numpy as np
import pandas as pd
from math import log, exp
import statsmodels.api as sm
import scipy
import sys

fileout,prefix,suffix,st,en,slide,columnname,covariance_length,chrlen,pval=sys.argv[1:]

covariance_length = int(float(covariance_length))
stepsize = float(slide)

plotting_covs = []

final_goodbyes = []
final_covs = []

keep_all_data = []

mns = []
sds = []


ctr = 0

for chromosome in range(int(float(st)),int(float(en))+1):
    try:
        table = pd.read_csv(prefix+str(chromosome)+suffix,sep='\t')
        subtable=table
        keep_all_data += subtable[columnname].to_list()
    except:
        pass
stddev = np.std(keep_all_data)
avg = np.mean(keep_all_data)


for chromosome in range(int(float(st)),int(float(en))+1):
    try:

        table = pd.read_csv(prefix+str(chromosome)+suffix,sep='\t')

        subtable=table
        #subtable['STANDARDIZED'] = ( subtable[columnname] - avg ) / stddev
        subtable['STANDARDIZED'] = ( subtable[columnname] - subtable[columnname].mean()) / subtable[columnname].std()

        mns.append(subtable[columnname].mean())
        sds.append(subtable[columnname].std())
        
        K = subtable.shape[0]
        xs = subtable['STANDARDIZED']

        byebye = covariance_length
        covs = []
        cors = []
        byes = []
        
        for bye in range(1, byebye + 1):
            seq1 = np.arange(1, K - bye + 1, bye)
            seq2 = np.arange(1 + bye, K + 1, bye)
            xs1 = xs[seq1-1]  # Adjusting index to 0-based for Python
            xs2 = xs[seq2-1]
            covs.append(np.cov(xs1, xs2)[0, 1])
            byes.append(bye)

        plotting_covs.append(covs)

        # plt.plot(byes,covs)

        final_goodbyes += byes
        final_covs += covs

        ctr += 1

    except:
        pass
    
final_goodbyes = np.array(final_goodbyes)
final_covs = np.array(final_covs)
final_goodbyes = final_goodbyes[final_covs>0.]
final_covs = final_covs[final_covs>0.]
    

xvar = np.array(final_goodbyes) * stepsize
# xvar = np.array(byes)
ytilde = - np.log(final_covs)

# xvar = sm.add_constant(xvar)  # Adds a constant term to the predictor
fittilde = sm.OLS(ytilde, xvar).fit()

print(fittilde.summary())

print('\n\n\n\n')

print('theta')
print(fittilde.params[0])

print('\n')
print('overall summary stats')
print(np.mean(keep_all_data))
print(np.std(keep_all_data))

print('\n')
print('meta analysis summary stats')
print(np.mean(mns))
print(np.mean(sds))

print('\n\n\n\n')

print(pd.Series(keep_all_data).describe())


### determining the upper bound in t-test ###

from scipy.stats import norm
from scipy.optimize import root_scalar

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

pval = float(pval)
siglevel = pval
chrnum = ctr
chrlen = float(chrlen)

theta = fittilde.params[0]
#stddev = np.std(keep_all_data)
#avg = np.mean(keep_all_data)
stddev = np.mean(sds)
avg = np.mean(mns)

# Using scipy's root_scalar to find the roots similar to R's uniroot
lin_root = root_scalar(OU_approx_cont, args=(theta, chrnum * chrlen, chrnum, siglevel), bracket=[1, 20])
lin = lin_root.root if lin_root.converged else None

lin2_root = root_scalar(OU_approx, args=(theta, stepsize, chrnum * chrlen, chrnum, siglevel), bracket=[1, 20])
lin2 = lin2_root.root if lin2_root.converged else None

if lin is not None:
    result_lin = lin * stddev + avg
    print(result_lin)
else:
    print("Root not found for OU_approx_cont")

if lin2 is not None:
    result_lin2 = lin2 * stddev + avg
    print(result_lin2)
else:
    print("Root not found for OU_approx")

bonferroni=pval/len(keep_all_data)
bz=norm.ppf(1-bonferroni)
bu=bz * stddev + avg

keep_all_data = np.array(keep_all_data)

def max_slide(x,step):
    J = len(x)-step+1
    y = []
    for j in range(J):
        mx = max(x[j:(j+step)])
        y.append(mx)
    return y

def min_slide(x,step):
    J = len(x)-step+1
    y = []
    for j in range(J):
        mn = min(x[j:(j+step)])
        y.append(mn)
    return y

f = open(fileout,'w')
f.write('p-value:'); f.write(str(pval)); f.write('\n')
f.write('num-tests:'); f.write(str(len(keep_all_data))); f.write('\n')
f.write('theta:'); f.write(str(theta)); f.write('\n')
f.write('mean:'); f.write(str(avg)); f.write('\n')
f.write('standard-deviation:'); f.write(str(stddev)); f.write('\n')
f.write('upper-discrete-z:'); f.write(str(lin2)); f.write('\n')
f.write('upper-continuous-z:'); f.write(str(lin)); f.write('\n')
f.write('upper-bonferroni-z:'); f.write(str(bz)); f.write('\n')
f.write('upper-discrete-raw:'); f.write(str(result_lin2)); f.write('\n')
f.write('upper-continuous-raw:'); f.write(str(result_lin)); f.write('\n')
f.write('upper-bonferroni-raw:'); f.write(str(bu)); f.write('\n')
f.write('upper-discrete-perc:'); f.write(str((keep_all_data>=result_lin2).mean())); f.write('\n')
f.write('upper-continuous-perc:'); f.write(str((keep_all_data>=result_lin).mean())); f.write('\n')
f.write('upper-bonferroni-perc:'); f.write(str((keep_all_data>=bu).mean())); f.write('\n')
f.write('upper-discrete-num-efftests-1:'); f.write(str((pd.Series(keep_all_data).rolling(window=1,step=1).max()>=result_lin2).sum())); f.write('\n')
f.write('upper-discrete-num-efftests-5:'); f.write(str((pd.Series(keep_all_data).rolling(window=5,step=5).max()>=result_lin2).sum())); f.write('\n')
f.write('upper-discrete-num-efftests-10:'); f.write(str((pd.Series(keep_all_data).rolling(window=10,step=10).max()>=result_lin2).sum())); f.write('\n')
f.write('upper-discrete-num-efftests-15:'); f.write(str((pd.Series(keep_all_data).rolling(window=15,step=15).max()>=result_lin2).sum())); f.write('\n')
f.write('upper-discrete-num-efftests-20:'); f.write(str((pd.Series(keep_all_data).rolling(window=20,step=20).max()>=result_lin2).sum())); f.write('\n')
f.write('upper-continuous-num-efftests-5:'); f.write(str((pd.Series(keep_all_data).rolling(window=5,step=5).max()>=result_lin).sum())); f.write('\n')
f.write('upper-continuous-num-efftests-10:'); f.write(str((pd.Series(keep_all_data).rolling(window=10,step=10).max()>=result_lin).sum())); f.write('\n')
f.write('upper-bonferroni-num-efftests-5:'); f.write(str((pd.Series(keep_all_data).rolling(window=5,step=5).max()>=bu).sum())); f.write('\n')
f.write('upper-bonferroni-num-efftests-10:'); f.write(str((pd.Series(keep_all_data).rolling(window=10,step=10).max()>=bu).sum())); f.write('\n')
f.write('actual-max:'); f.write(str(max(keep_all_data))); f.write('\n')
keep_all_data = list(keep_all_data)
f.write('actual-minmax-of-2:'); f.write(str(max(min_slide(keep_all_data,2)))); f.write('\n')
f.write('actual-minmax-of-3:'); f.write(str(max(min_slide(keep_all_data,3)))); f.write('\n')
f.write('actual-minmax-of-4:'); f.write(str(max(min_slide(keep_all_data,4)))); f.write('\n')
f.write('actual-minmax-of-5:'); f.write(str(max(min_slide(keep_all_data,5)))); f.write('\n')
f.close()
