# RCAnalysis

Work in progress

## RCAnalysis_SplitReads
Takes as input MinION fastq files and primer sequence
Searches for the primer sequence in each read cuts them up, creating a fastq file for each
  Has option for RCA and dumbell methods
TODO: set filter for read length after cutting

## consensus.sh
Takes as input the cut up fastq files for each read and creates a consensus for each
Two methods available
Output is a fastq (or fasta) file with the consensus reads (one for each *amplicon*) that can be used downstream as usual