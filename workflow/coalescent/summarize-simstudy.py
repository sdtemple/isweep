#!/bin/python

# require a python version, environment with pandas, numpy

# Seth D. Temple, sdtemple@uw.edu
# April 19, 2024
# Simple scripting to summarize simulation study on coalescent IBD

# Example terminal command:
# python summarize-simstudy.py filein fileout column
# filein : a tab-delimited file
# fileout : prefix for fileout.tsv
# column : integer; start with 1; script will convert to python indexing
# typ : string; type of simulation study; give it a name

# imports
import sys
import pandas as pd
import numpy as np

# inputs
filein, fileout, column, typ = sys.argv[1:]
column=int(float(column))-1
tab = pd.read_csv(filein, sep='\t', header=None)

# setup
uperturb = tab[column].unique().tolist()
fixedleft = tab.loc[0,5:(column-1)]
fixedright = tab.loc[0,(column+1):]
fixed = fixedleft.tolist() + fixedright.tolist()
truth = tab[0].unique().tolist()

# file to write to
g=open(fileout, 'w')
g.write('TYPE\tPARAM\tTRUE\tAVG\tBIAS\tLOW\tAVG\tUPP\tCOV\n') # modify if you add more statistics

# loop over parameter perturbations
for tru in truth:
    for perturb in uperturb:
        subleft = tab[tab[column]==perturb]

        # here define statistics of interest
        # column 0 is true selection coefficient
        # column 1 is inital estimate selection coefficient
        # column 2 is bootstrap (bias-corrected) median selection coefficient
        # column 3 is bootstrap (bias-corrected) left bound selection coefficient
        # column 4 is bootstrap (bias-corrected) right bound selection coefficient

        avg = subleft[2].mean() # mean over replicates of bootstrap median
        low = subleft[3].mean() # mean over replicates of bootstrap left bound
        upp = subleft[4].mean() # mean over replicates of bootstrap right bound
        bia = tru - avg # mean bias
        cov = ((subleft[3]<=tru)&(subleft[4]>=tru)).mean() # coverage

        g.write(typ+'\t')
        g.write(str(perturb)+'\t')
        g.write(str(tru)+'\t')
        g.write(str(avg)+'\t')
        g.write(str(bia)+'\t')
        g.write(str(low)+'\t')
        g.write(str(avg)+'\t') # duplicate avg
        g.write(str(upp)+'\t')
        g.write(str(cov)+'\n') # make sure final one has \n new line

g.close()

# ### For selection coefficient simulation study ###
#
# # loop over parameter perturbations
# for tru in truth:
#     subleft = tab[tab[column]==tru]
#
#     # here define statistics of interest
#     # column 0 is true selection coefficient
#     # column 1 is inital estimate selection coefficient
#     # column 2 is bootstrap (bias-corrected) median selection coefficient
#     # column 3 is bootstrap (bias-corrected) left bound selection coefficient
#     # column 4 is bootstrap (bias-corrected) right bound selection coefficient
#
#     avg = subleft[2].mean() # mean over replicates of bootstrap median
#     low = subleft[3].mean() # mean over replicates of bootstrap left bound
#     upp = subleft[4].mean() # mean over replicates of bootstrap right bound
#     bia = tru - avg # mean bias
#     cov = ((subleft[3]<=tru)&(subleft[4]>=tru)).mean() # coverage
#
#     g.write(typ+'\t')
#     g.write(str(tru)+'\t')
#     g.write(str(tru)+'\t')
#     g.write(str(avg)+'\t')
#     g.write(str(bia)+'\t')
#     g.write(str(low)+'\t')
#     g.write(str(avg)+'\t') # duplicate avg
#     g.write(str(upp)+'\t')
#     g.write(str(cov)+'\n') # make sure final one has \n new line
#
# g.close()
