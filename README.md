# RCAnalysis

Work in progress

## RCAnalysis_SplitReads
Takes as input MinION fastq files and primer sequence
Searches for the primer sequence in each read cuts them up, creating a fastq file for each
  Has option for RCA and dumbell methods
  It does not remove the primer sequence, so these must be sof clipped in alignment (or trimmed before)
    Does this cause issue with the consensus calling? I dont think so (but its easy to edit if it does)


## consensus.sh
Takes as input the cut up fastq files for each read and creates a consensus for each
  make sure only taking those with multiple reads, and then adding the singletons at the end if user wants
Two methods available
Output is a fastq (or fasta) file with the consensus reads (one for each *amplicon*) that can be used downstream as usual