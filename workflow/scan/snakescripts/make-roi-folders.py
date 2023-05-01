#!/bin/python

# Make folders for regions of interest
# Seth D. Temple, sdtemple@uw.edu
# April 28, 2023

# setup
import pandas as pd

roitable=snakemake.input.roitable
study=snakemake.config['CHANGE']['FOLDERS']['STUDY']
micro=roitable
macro=study

sims = pd.read_csv(micro, sep='\t')
# sims = pd.read_csv(micro, sep='\t', header=0)
J = sims.shape[0]
for j in range(J):
	row = sims.loc[j,]
	if not os.path.exists(macro+'/roi'+str(row.NAME)):
		os.mkdir(macro+'/roi'+str(row.NAME))
	if not os.path.exists(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)):
		os.mkdir(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM))
	if not os.path.exists(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)):
		os.mkdir(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER))
	if not os.path.exists(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)+ '/left'+str(row.BPLEFT)):
		os.mkdir(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)+ '/left'+str(row.BPLEFT))
	if not os.path.exists(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)+ '/left'+str(row.BPLEFT)+ '/right'+str(row.BPRIGHT)):
		os.mkdir(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+ '/center'+str(row.BPCENTER)+ '/left'+str(row.BPLEFT)+ '/right'+str(row.BPRIGHT))
	if not os.path.exists(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)+ '/left'+str(row.BPLEFT)+ '/right'+str(row.BPRIGHT)+'/ready'):
		os.mkdir(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+ '/center'+str(row.BPCENTER)+ '/left'+str(row.BPLEFT)+ '/right'+str(row.BPRIGHT)+'/ready')

h.open(study+'/roifolders.done','w')
h.write('Hello world!')
h.close()
