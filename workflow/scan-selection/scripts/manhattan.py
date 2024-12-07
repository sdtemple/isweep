# make a manhattan plot of IBD rates

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# input parameters
filein = sys.argv[1]
fileout = sys.argv[2]
chrlow = int(sys.argv[3])
chrhigh = int(sys.argv[4])
cutoff1 = float(sys.argv[5])
cutoff2 = float(sys.argv[6])
title = sys.argv[7]
height = sys.argv[9]
width = sys.argv[8]
height = float(height)
width = float(width)

# reading in data
ibd = pd.read_table(filein, sep='\t')

# maths
medi = ibd["COUNT"].median()
stdv = ibd["COUNT"].std()
a = medi - stdv * cutoff2
b = medi + stdv * cutoff2
low = a
newibd = ibd[(ibd["COUNT"] >= a) & (ibd["COUNT"] <= b)].copy()

# compute browning and browning 2020 cutoff
medi = newibd["COUNT"].median(numeric_only=True)
stdv = newibd["COUNT"].std(numeric_only=True)
upp = medi + stdv * cutoff1
md = medi

# chromosome maths
ibd["CUMPOS"] = None
s = 0
mbp = [0] * (chrhigh + 1)
chr = {i:i for i in range(chrlow, chrhigh + 1)}
nchr = len(chr)

# ibd=newibd
for i in range(chrlow, chrhigh + 1):
    mbp[i] = ibd.loc[ibd["CHROM"] == chr[i], "CMWINDOW"].max()
    ibd.loc[ibd["CHROM"] == chr[i], "CUMPOS"] = ibd.loc[ibd["CHROM"] == chr[i], "CMWINDOW"] + s
    s = s + mbp[i]
ibd["CHROM"] = ibd["CHROM"].astype(int)

# plotting
plt.figure(figsize=(width, height))
mxy=ibd['COUNT'].max()*1.1
mny=ibd['COUNT'].min()
plt.ylim(mny,mxy)
for chrom, group in ibd.groupby("CHROM"):
    if chrom % 2 == 0:
        clr = "black"
    if chrom % 2 == 1:
        clr = "gray"
    # if dots are preferred
    # plt.scatter(group["CUMPOS"], group["COUNT"], c=clr, s=0.1)
    # if lines are preferred
    plt.plot(group["CUMPOS"], group["COUNT"], c=clr)

plt.axhline(y=md, color="blue", linestyle="--", label="median")
plt.axhline(y=upp, color="red", linestyle="--", label='cutoff')
plt.ylabel("IBD rate")
plt.yticks([])
plt.ylim(mny, mxy)
plt.xlim(ibd['CUMPOS'].min()-25,ibd['CUMPOS'].max()+25)
plt.title(title, loc='left')
plt.tight_layout()

plt.legend(loc="upper right",title=None)

# Customize x-axis ticks
# Get unique chromosome values
chromosome_values = sorted(ibd["CHROM"].unique())
x_ticks = [np.mean(ibd[ibd["CHROM"] == chrom]["CUMPOS"]) for chrom in chromosome_values]
x_labels = [f"{chrom}" for chrom in chromosome_values]
plt.xticks(x_ticks, x_labels, rotation=90)

# saving plots
pic = "jpeg"
plt.savefig(f"{fileout}.{pic}")
pic = "png"
plt.savefig(f"{fileout}.{pic}")
# pic='tiff'
# plt.savefig(f"{fileout}.{pic}")
# pic='eps'
# plt.savefig(f"{fileout}.{pic}")
