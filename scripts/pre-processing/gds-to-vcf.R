library(SeqArray)
library(SNPRelate)
library(gdsfmt)
args=commandArgs(trailingOnly=TRUE)

# Via SeqArray
filein=seqOpen(args[1])
fileout=args[2]
seqGDS2VCF(filein, fileout, info.var=character(0))
seqClose(filein)

# # via SeqArray
# tryCatch({
#   filein=seqOpen(args[1])
#   fileout=args[2]
#   seqGDS2VCF(filein, fileout, info.var=character(0))
#   seqClose(args[1])
# }, error=function(e) {
#   print(e)
#   seqSNP2GDS(args[1],paste(args[1],".seq.gds",sep=""))
#   filein=seqOpen(paste(args[1],".seq.gds",sep=""))
#   fileout=args[2]
#   seqGDS2VCF(filein, fileout, info.var=character(0))
#   seqClose(filein)
# })