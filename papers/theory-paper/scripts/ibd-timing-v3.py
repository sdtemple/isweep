from isweep2 import *
import time
import sys
import networkx as nx

fileout, Ne_file, nboot, cM = sys.argv[1:]

f = open(fileout,'w')
f.write('alg\t2k\t4k\t8k\t16k\t32k\t64k\t128k\n')

cM = float(cM)

afr_Ne = read_Ne(Ne_file)
num = int(float(nboot))
nboot = float(nboot)

# full algorithm
f.write('both')
f.write('\n')

ns = [2_000,4_000,8_000,16_000,32_000,64_000,128_000]
# ns = [2_000,4_000]

timeds = []
for i in range(len(ns)):
    n = ns[i]
    print(n)
    timed = []
    for j in range(num):
        print(j)
        start = time.time()
        out = simulate_ibd_timing(n,
                                  afr_Ne,
                                  long_ibd=cM,
                                  short_ibd=cM,
                                  pairwise_output=False,
                                  record_dist=False,
                                  execute_merge=True,
                                  execute_prune=True,
                                 )
        end = time.time()
        delta = end - start
        timed.append(delta)
    print('')
    timeds.append(timed)
    
for ts in timeds[:-1]:
    sm = sum(ts) / nboot
    f.write(str(sm))
    f.write('\t')
sm = sum(timeds[-1]) / nboot
f.write(str(sm))
f.write('\n')

# w/o prune

f.write('woprune\t')

ns = [2_000,4_000,8_000,16_000,32_000,64_000,128_000]
# ns = [2_000,4_000]

timeds = []
for i in range(len(ns)):
    n = ns[i]
    print(n)
    timed = []
    for j in range(num):
        print(j)
        start = time.time()
        out = simulate_ibd_timing(n,
                                  afr_Ne,
                                  long_ibd=cM,
                                  short_ibd=cM,
                                  pairwise_output=False,
                                  record_dist=False,
                                  execute_merge=True,
                                  execute_prune=False
                                 )
        end = time.time()
        delta = end - start
        timed.append(delta)
    print('')
    timeds.append(timed)
    
for ts in timeds[:-1]:
    sm = sum(ts) / nboot
    f.write(str(sm))
    f.write('\t')
sm = sum(timeds[-1]) / nboot
f.write(str(sm))
f.write('\n')

# w/o merge

f.write('womerge\t')

ns = [2_000,4_000,8_000,16_000,32_000,64_000,128_000]
# ns = [2_000,4_000]

timeds = []
for i in range(len(ns)):
    n = ns[i]
    print(n)
    timed = []
    for j in range(num):
        print(j)
        start = time.time()
        out = simulate_ibd_timing(n,
                                  afr_Ne,
                                  long_ibd=cM,
                                  short_ibd=cM,
                                  pairwise_output=False,
                                  record_dist=False,
                                  execute_merge=False,
                                  execute_prune=True
                                 )
        end = time.time()
        delta = end - start
        timed.append(delta)
    print('')
    timeds.append(timed)
    
for ts in timeds[:-1]:
    sm = sum(ts) / nboot
    f.write(str(sm))
    f.write('\t')
sm = sum(timeds[-1]) / nboot
f.write(str(sm))
f.write('\n')

# w/o both

f.write('woboth\t')

ns = [2_000,4_000,8_000]
# ns = [2_000,4_000]

timeds = []
for i in range(len(ns)):
    n = ns[i]
    print(n)
    timed = []
    for j in range(num):
        print(j)
        start = time.time()
        out = simulate_ibd_timing(n,
                                  afr_Ne,
                                  long_ibd=cM,
                                  short_ibd=cM,
                                  pairwise_output=False,
                                  record_dist=False,
                                  execute_merge=False,
                                  execute_prune=False
                                 )
        end = time.time()
        delta = end - start
        timed.append(delta)
    print('')
    timeds.append(timed)
    
for ts in timeds[:-1]:
    sm = sum(ts) / nboot
    f.write(str(sm))
    f.write('\t')
sm = sum(timeds[-1]) / nboot
f.write(str(sm))
f.write('\n')

# independent

f.write('indep\t')

ns = [2_000,4_000,8_000,16_000,32_000]
# ns = [2_000,4_000]

timeds = []
for i in range(len(ns)):
    n = ns[i]
    print(n)
    timed = []
    for j in range(num):
        print(j)
        start = time.time()
        out = simulate_ibd_independent(n,
                                       afr_Ne,
                                       long_ibd=cM
                                      )
        end = time.time()
        delta = end - start
        timed.append(delta)
    print('')
    timeds.append(timed)
    
for ts in timeds[:-1]:
    sm = sum(ts) / nboot
    f.write(str(sm))
    f.write('\t')
sm = sum(timeds[-1]) / nboot
f.write(str(sm))
f.write('\n')

# erdos renyi

f.write('erdosrenyi\t')

ns = [2_000,4_000,8_000,16_000,32_000]
# ns = [2_000,4_000]

timeds = []
for i in range(len(ns)):
    n = ns[i]
    print(n)
    timed = []
    sample_size2 = int(n * 2)
    for j in range(num):
        print(j)
        start = time.time()
        prop = probability_ibd_isweep(0.0, 1., afr_Ne, cM)
        the_graph = nx.erdos_renyi_graph(sample_size2, prop)
        end = time.time()
        delta = end - start
        timed.append(delta)
    print('')
    timeds.append(timed)
    
for ts in timeds[:-1]:
    sm = sum(ts) / nboot
    f.write(str(sm))
    f.write('\t')
sm = sum(timeds[-1]) / nboot
f.write(str(sm))
f.write('\n')
