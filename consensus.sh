#!/bin/bash


# required software
# R and packages
# racon
# medaka
# minimap2
# filtlong


method = "A"
file_dir_in = "path to reads"

# Only method B options
RNDracon = 3
threads = 1
mod = "r941_min_fast_g330"
# this takes in the deconcatenated fastq files and calls a consensus that can then be used downstream for alignment
# below is based on my MinION_UMI script

## Method A ##
if [[ $method == "A" || $method == "a" || $method == "C" || $method == "c" ]]
then
    Rscript "${script_dir}/methodA.R.R" "$file_dir_in"
    # bowtie2 --local -x ${ref/.fa*/} -U "${file_in/.f*/.consensus.fastq}" -S "${name}.methodB.consensus.sam" --threads $threads
    # samtools view -bS "${name}.methodB.consensus.sam" | samtools sort -o "${name}.methodB.consensus.bam"
    # samtools index "${name}.methodB.consensus.bam"
fi


if [[ $method == "B" || $method == "b" || $method == "C" || $method == "c" ]]
then
    ##############
    ## Method B ##

    # this looses the quality information, so polish before?

    # if only intresed in SNPs and not indels, then use the actual reference,
    # otherwise, this makes more sense to me? rethinkg this
    # you could take the "best read" from the umis and use that as the initial reference?

    count=1
    for reads in $(ls $file_dir_in/*.fq)
    do
        #this is a reference agnostic approach
        filtlong --keep_percent 1 "$reads" | head -n 2 | tail -n 1 >> "read_${count}_best.fa"
        i=0
        while [[ $i -lt $RNDracon ]]
        do
            draft="read_${count}_best.fa" #to readsm above - will get overwritten here
            minimap2 -ax map-ont "${draft}" "$reads" > "read_${count}.sam"
            aln="read_${count}.sam"

            if [[ $racon_or_medaka == "racon" ]]
            then
                racon -w 500 -m 8 --quality-threshold 10 --error-threshold 0.3 -x -1 -g -4 -q -1 -t ${threads} \
                  "${reads}" "${aln}" "${draft}" > "${draft}.tmp"
                mv ${draft}.tmp ${draft}
            fi

            # Medaka has been trained to correct draft sequences processed through racon, specifically racon run four times iteratively with: racon -m 8 -x -6 -g -8 -w 500 ...
            if [[ $i == $RNDracon || $racon_or_medaka != "racon" ]]
            then
                medaka_consensus --model $mod -i "${reads}" -d "${draft}" -t ${threads} -o $( dirname $draft )
                mv "$( dirname $draft )/consensus.fasta" ${draft}
            fi
            i=$[$i+1]
        done
    done

    cat ${file_dir_in}/*_best.fa >> "methodB.consensus.fa"
    
    # bowtie2 --local -x ${ref/.fa*/} -f "${name}.methodA.consensus.fa" -S "${name}.methodA.consensus.sam" --threads $threads
    # samtools view -bS "${name}.methodA.consensus.sam" | samtools sort -o "${name}.methodA.consensus.bam"
    # samtools index "${name}.methodA.consensus.bam"
fi
