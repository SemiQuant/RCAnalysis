require(muscle)
require(Biostrings)
require(tidyverse)
require(sarlacc)

args = commandArgs(trailingOnly=TRUE)
dir_in <- args[1]
# dir_in <- "/Users/SemiQuant/Bioinformatics/Projects/RCAnalysis/test"
files <- list.files(path = dir_in, pattern = "_multiple.fq", full.names = T)


for (reads_in in files){
  read.seq <- readQualityScaledDNAStringSet(reads_in)
  msa.out <- multiReadAlign(read.seq, rep(1, length(read.seq)))
  # We create consensus sequences from these MSAs, representing the error-corrected sequence.
  # The quality scores are constructed from the qualities of the individual read sequences.
  # Higher consensus qualities for a position indicate that many reads are in agreement.
  cons.out <- consensusReadSeq(msa.out)
  writeXStringSet(cons.out, paste0(dir_in, "methodA_consensus.fastq"),
                  append=T, format="fastq", compress = F, qualities = cons.out@quality)
}
