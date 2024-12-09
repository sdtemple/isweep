# Make a plot for modeling the selective sweep.
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-07-22
# Description: Plot the mean and quantiles of allele frequency trajectories under a selective sweep

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from math import *
from scipy.stats import binom
from scipy.stats import norm
import argparse

# Argument parser setup
parser = argparse.ArgumentParser(description='Plot a selective sweep with uncertainty')
parser.add_argument('--output_file', 
                    type=str,
                    required=True, 
                    help='Output file for the plot')
parser.add_argument('--s', 
                    type=float, 
                    required=True,
                    help='Selection coefficient')
parser.add_argument('--su', 
                    type=float,
                    required=True, 
                    help='Selection coefficient upper limit')
parser.add_argument('--z', 
                    type=float,
                    required=True, 
                    help='Z value for normal distribution -based confidence intervals')
parser.add_argument('--p', 
                    type=float,
                    required=True, 
                    help='Present-day allele frequency')
parser.add_argument('--Ne', 
                    type=str,
                    required=True, 
                    help='Effective population sizes file')
parser.add_argument('--standing_variation', 
                    type=float, 
                    default=-0.01, 
                    help='(default: -0.01) Standing variation frequency')
parser.add_argument('--genetic_model', 
                    type=str, 
                    default='a', 
                    help='(default: additive) Genetic selection model')
parser.add_argument('--ploidy', 
                    type=int, 
                    default=2, 
                    help='(default: 2) Ploidy level')
parser.add_argument('--nboot', 
                    type=int, 
                    default=1000, 
                    help='(default: 1000) Number of Wright-Fisher bootstraps')
parser.add_argument('--upper_quantile', 
                    type=float, 
                    default=0.99, 
                    help='(default: 0.99) Upper quantile of bootstraps')
parser.add_argument('--lower_quantile', 
                    type=float, 
                    default=0.01, 
                    help='(default: 0.01) Lower quantile of bootstraps')
parser.add_argument('--xaxis_length', 
                    type=int, 
                    default=150, 
                    help='(default: 150) Length of the time axis (generations)')
parser.add_argument('--line_color', 
                    type=str, 
                    default='k', 
                    help='(default: black) Color for plotting')
parser.add_argument('--font_size', 
                    type=int, 
                    default=14, 
                    help='(default: 14) Font size for plotting')
parser.add_argument('--alpha', 
                    type=float, 
                    default=0.33, 
                    help='Alpha for grid lines in plotting')
parser.add_argument('--title', 
                    type=str, 
                    default=None, 
                    help='Title for the plot')

args = parser.parse_args()

plt.rc("font", size=args.font_size)

nboot = args.nboot
length = args.xaxis_length
xax = [i for i in range(length)]

s = args.s
p = args.p
su = args.su
the_color = args.line_color

def mean(lst):
    return sum(lst) / len(lst)

def walk_variant_backward(s, p0, Ne, random_walk=False, one_step_model='a', tau0=0, sv=-0.01, ploidy=2):
    assert ploidy in [1, 2]
    assert one_step_model in ['m', 'a', 'r', 'd']
    assert p0 <= 1
    assert p0 >= 0
    assert sv < 1

    def haploid_bwd(p, s): 
        return p / (1 + s - s * p)

    def multiplicative_bwd(p, s):
        num = 1
        dnm = 1 + (1 - p) * s
        return p * num / dnm

    def dominant_bwd(p, s):
        a = p * s
        b = 1 + s - 2 * p * s
        c = -p
        qf = -b + sqrt((b ** 2) - 4 * a * c)
        qf = qf / 2 / a
        return qf

    def additive_bwd(p, s):
        a = s
        b = 1 + s - 2 * p * s
        c = -p
        qf = -b + sqrt((b ** 2) - 4 * a * c)
        qf = qf / 2 / a
        return qf

    def recessive_bwd(p, s):
        a = (1 - p) * s
        b = 1
        c = -p
        qf = -b + sqrt((b ** 2) - 4 * a * c)
        qf = qf / 2 / a
        return qf

    if ploidy == 1:
        one_step = haploid_bwd
    else:
        if one_step_model == 'a':
            one_step = additive_bwd
        elif one_step_model == 'r':
            one_step = recessive_bwd
        elif one_step_model == 'd':
            one_step = dominant_bwd
        else:
            one_step = multiplicative_bwd

    ps = []  
    xs = []  
    Ns = []  
    t = floor(tau0)
    p = p0
    N = Ne[0]
    x = floor(p * ploidy * N)
    Ns.append(N)
    xs.append(x)
    ps.append(p)

    if random_walk:  
        for G in range(1, max(Ne.keys()) + 1):
            try: 
                N = Ne[G]
            except KeyError:
                pass
            if G > t:
                p = one_step(p, s)
            x = int(binom.rvs(int(ploidy * N), p))
            p = x / ploidy / N
            if x < 1:
                break
            if p >= 1:
                break
            if p <= sv:
                s = 0
            ps.append(p)
            xs.append(x)
            Ns.append(N)

        return np.array([ps, Ns, xs], dtype=float)  

    else:  
        for G in range(1, max(Ne.keys()) + 1):
            try: 
                N = Ne[G]
            except KeyError:
                pass
            if G > t:
                p = one_step(p, s)
            x = floor(p * ploidy * N)
            if x < 1:
                break
            if p >= 1:
                break
            if p <= sv:
                s = 0
            Ns.append(N)
            xs.append(x)
            ps.append(p)

        return np.array([np.array(ps), np.array(Ns), np.array(xs)])  

def read_Ne(file):
    Ne = {}
    with open(file, 'r') as f:
        f.readline() 
        for line in f:
            g, size = line.split('\t')[:2]
            Ne[int(g)] = int(float(size))
    return Ne

Ne = read_Ne(args.Ne)

yax = [[] for y in range(len(xax))]

a, b, c = walk_variant_backward(s, p, Ne, one_step_model=args.genetic_model,ploidy=args.ploidy,sv=args.standing_variation)
A = a[:length]

sig = (su - s) / args.z

for b in range(nboot):
    try:
        sb = norm.rvs(s, sig)
        a, b, c = walk_variant_backward(sb, p, Ne, random_walk=True, one_step_model=args.genetic_model,ploidy=args.ploidy,sv=args.standing_variation)
        A = a[:length]
        for j in range(len(xax)):
            yax[j].append(A[j])
    except:
        pass

A = [mean(y) for y in yax]
plt.plot(xax, A, color=the_color, label='Mean', linewidth=2)

A = [np.quantile(y, args.upper_quantile) for y in yax]
plt.plot(xax, A, color=the_color, label='Upper percentile', linewidth=2, linestyle='dashed')

A = [np.quantile(y, args.lower_quantile) for y in yax]
plt.plot(xax, A, color=the_color, label='Lower percentile', linewidth=2, linestyle='dotted')

plt.legend(title=None, loc="upper right", shadow=True)
plt.ylim(0, 1)
plt.xlabel("Time (generations)")
plt.ylabel('Allele frequency')
plt.title(args.title)
plt.grid(alpha=args.alpha)
plt.savefig(args.output_file)
