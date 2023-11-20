#!/usr/bin/env Rscript
library(gdsfmt)
library(SNPRelate)
args=commandArgs(trailingOnly=TRUE)
snpgdsVCF2GDS(vcf.fn=args[1],out.fn=args[2],method=args[4],ignore.chr.prefix=args[3])
