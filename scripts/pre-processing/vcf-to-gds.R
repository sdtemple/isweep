library(SeqArray)
args=commandArgs(trailingOnly=TRUE)

# via SeqArray
filein=args[1]
fileout=args[2]
info.import=c("AC","AF","DP")
seqVCF2GDS(filein, fileout, info.import=info.import)