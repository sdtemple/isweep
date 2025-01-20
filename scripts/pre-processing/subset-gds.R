library(SeqArray)

### command line arguments

args <- commandArgs(trailingOnly=T)
input.file <- args[1]
subset.file <- args[2]
mac <- as.integer(args[3])
missingness <- as.numeric(args[4])
outputfileprefix <- args[5]

outgds <- paste(outputfileprefix,".gds",sep="")
outvcf <- paste(outputfileprefix,".vcf.gz",sep="")
outsum <- paste(outputfileprefix,".summary.stats.tsv",sep="")
# outvar <- paste(outputfileprefix,".kept.variants.txt",sep="")
# outsum <- paste(outputfileprefix,".summary.stats.rds",sep="")
# outvar <- paste(outputfileprefix,".kept.variants.rds",sep="")

### summary stats, courtesy liz blue ###

sumstat <- data.frame(variant.id  = seqGetData(input.file, "variant.id"), mac = seqAlleleCount(input.file, ref.allele=0L, minor=T), maf = seqAlleleFreq(input.file, ref.allele=0L, minor=T), miss = seqMissing(input.file, TRUE))
my.variants <- subset(sumstat$variant.id, sumstat$mac >= mac, sumstat$miss < missingness)
write.table(sumstat, file=outsum, sep="\t", row.names=F)
# writeLines(my.variants, con=outvar)
# saveRDS(my.variants, file=outvar)
# saveRDS(sumstat, file=outsum)

### subset gds file ###

gdsfile <- seqOpen(input.file)
print("opened the gds file")

my.samples <- readLines(subset.file)
seqSetFilter(gdsfile, variant.id = my.variants, sample.id = my.samples)
print("subset the gds file by sample and variant")

### export to gds and/or vcf ###

seqGDS2VCF(gdsfile,outvcf,info.var=character(0))
print("saved the filtered vcf file")

# seqExport(gdsfile, outgds)
# print("saved the filtered gds file")

seqClose(gdsfile)

### the end ###

