import sys
from isweep import *

fileout = sys.argv[1]
sample_size = sys.argv[2]
constant_Ne = sys.argv[3]
ibd_length = sys.argv[4]
num_sim = sys.argv[5]

f = open(fileout,'w')

N = int(float(constant_Ne))
n = int(float(sample_size))
W = float(ibd_length)
w =  W / 100

ploidy = 2

q = (2 * ploidy * N * w / 2 + 1 ) ** (-1)
# divide by 2 b/c simulating from erlang
# times by two is from calculation
m = ploidy * n
mq = m * (m-1) / 2 * q

f.write('N=' + str(N) + '\n')
f.write('n=' + str(n) + '\n')
f.write('w=' + str(w) + '\n')
f.write('E=' + str(mq) + '\n')
f.write('q=' + str(q) + '\n')
f.write('-----\n')
f.write('num_tracts\tif_singleton\n')

K = int(float(num_sim))
for k in range(K):
    data = simulate_ibd_constant(n, N, W, W, ploidy, True, True)
    nt = data[0]
    ns = sum([1 if (item[1]+item[0]==1) else 0 for item in data[-1]])
    f.write(str(nt) + '\t' + str(ns) + '\n')
f.close()
