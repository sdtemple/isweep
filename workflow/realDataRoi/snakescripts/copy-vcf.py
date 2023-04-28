# isweep real data analysis part before 2 (copy vcf from scan to roi)
# Seth D. Temple, sdtemple@uw.edu
# April 27, 2023

# making folder structure
import os
import sys
import pandas as pd
study,roi=sys.argv[1:]
# study is the name of a folder for a real data analysis
# roi is the name of a *.tsv file with info on an roi (region of interest)
# the header should be NAME CHROM BPLEFT BPCENTER BPRIGHT
# you can put cM position as well if you want to
# you do this after a sweep scan
macro=snakemake.config['CHANGE']['FOLDERS']['STUDY']
micro=snakemake.config['CHANGE']['FOLDERS']['ROI']
sims = pd.read_csv(micro, sep='\t')
sims.to_csv(macro+'/ROI.tsv',sep='\t')
J = sims.shape[0]
for j in range(J):
    row = sims.loc[j,]
    if not os.path.exists(macro+'/roi'+str(row.NAME)):
        os.mkdir(macro+'/roi'+str(row.NAME))
    if not os.path.exists(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)):
        os.mkdir(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM))
    if not os.path.exists(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)):
        os.mkdir(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER))
    if not os.path.exists(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)+'/left'+str(row.BPLEFT)):
        os.mkdir(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)+'/left'+str(row.BPLEFT))
    if not os.path.exists(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)+'/left'+str(row.BPLEFT)+'/right'+str(row.BPRIGHT)):
        os.mkdir(macro+'/roi'+str(row.NAME)+'/chr'+str(row.CHROM)+'/center'+str(row.BPCENTER)+'/left'+str(row.BPLEFT)+'/right'+str(row.BPRIGHT))
sims['FOLDER'] = [macro +'/roi'+str(sims.loc[j].NAME)+'/chr'+str(row.CHROM)+'/center'+str(sims.loc[j].BPCENTER)+'/left'+str(sims.loc[j].BPLEFT)+'/right'+str(sims.loc[j].BPRIGHT) for j in range(J)]

# copying vcfs
import shutil
chrlist=sims['CHROM'].tolist()
folderlist=sims['FOLDER'].tolist()
for i in range(len(chrlist)):
    chr=chrlist[i]
    folder=folderlist[i]
    print(folder)
    oldvcf=macro+'/vcfs/chr'+str(chr)+'.vcf.gz'
    newvcf=folder+'/chr'+str(chr)+'.vcf.gz'
    shutil.copy(oldvcf,newvcf)
