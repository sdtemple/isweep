library(SeqArray)

args <- commandArgs(trailingOnly=T)

gdsfile <- args[1]

number_to_keep <- as.integer(args[2])

outputfile <- args[3]

sample_ids <- seqGetData(gdsfile, "sample.id")

keep <- min(number_to_keep,length(sample_ids))

subset_of_ids <- sample_ids[1:keep]

writeLines(subset_of_ids, con=outputfile)
