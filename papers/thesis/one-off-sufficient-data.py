
import sys
from isweep import *

# adjust as desired
# bins = [2.0 + i * 0.1 for i in range(10)]
bins = [3.0 + i * 0.1 for i in range(10)]
bins += [4.0 + i * 0.2 for i in range(5)]
bins += [5.,5.5,6,np.inf]

fileout = sys.argv[1]
nboot = int(sys.argv[2])
Ne = read_Ne(sys.argv[4])
n = int(sys.argv[3])
p = float(sys.argv[5])
s = float(sys.argv[6])
m = n * 2
N = m * (m-1) / 2

f = open(fileout,'w')
f.write('selcoef\tlengthdistr\tlengthnum\tlabeledlengthdistr\tlabeledlengthnum\n')
for b in range(nboot):

    f.write(str(s)); f.write('\t')

    out = simulate_ibd_isweep(n,s,p,Ne,bins[0],bins[0])
    whole_length_distribution = big_format_distribution(out[0][2],out[0][4])
    whole_bins = bin_ibd_segments(whole_length_distribution,bins)

    whole_length_distribution_1 = big_format_distribution(out[1][2],out[1][4])
    whole_bins_1 = bin_ibd_segments(whole_length_distribution_1,bins)

    whole_length_distribution_0 = big_format_distribution(out[2][2],out[2][4])
    whole_bins_0 = bin_ibd_segments(whole_length_distribution_0,bins)

    val = minimize_scalar(chi2_isweep,
                args=(p,Ne,N,whole_bins,bins),
                bounds=(0,0.5),
                method='bounded'
               ).x
    f.write(str(val)); f.write('\t')

    val = minimize_scalar(chi2_isweep,
                args=(p,Ne,N,(sum(whole_bins),),(bins[0],np.inf)),
                bounds=(0,0.5),
                method='bounded'
               ).x
    f.write(str(val)); f.write('\t')

    val = minimize_scalar(chi2_labeled_isweep,
                args=(p,Ne,N,whole_bins_1,whole_bins_0,bins),
                bounds=(0,0.5),
                method='bounded'
               ).x
    f.write(str(val)); f.write('\t')

    val = minimize_scalar(chi2_labeled_isweep,
                args=(p,Ne,N,(sum(whole_bins_1),),(sum(whole_bins_0),),(bins[0],np.inf)),
                bounds=(0,0.5),
                method='bounded'
               ).x
    f.write(str(val)); f.write('\n')

f.close()

