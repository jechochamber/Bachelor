.libPaths("~/Rlibs")

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(Biostrings)


gtf <- '/fast/AG_Ohler/jdemoli/bachelorgit/data/Saccharomyces_cerevisiae.R64-1-1.113.gtf'
fafile <- '/fast/AG_Ohler/jdemoli/bachelorgit/data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'
fputrfile <- '/fast/AG_Ohler/jdemoli/bachelorgit/data/gencode_ccds_fputrs_human.fa'

stopifnot(file.exists(fafile))
stopifnot(file.exists(gtf))
fafileob = Rsamtools::FaFile(fafile)
Rsamtools::indexFa(fafile)

#get utrs
#load the gtf as a granges
gtf_gr <- rtracklayer::import(con=gtf,format='gtf')

#make a transcript db object
gtftxdb <- makeTxDbFromGRanges(gtf_gr)

#get utrs from this as grangeslists
fputrs <- fiveUTRsByTranscript(gtftxdb, use.names=TRUE)

#extract sequence and print to a file
GenomicFeatures::extractTranscriptSeqs(fputrs,x=fafileob) %>% writeXStringSet(fputrfile)



