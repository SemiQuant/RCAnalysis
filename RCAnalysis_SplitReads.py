#!/usr/bin/env python
# coding: utf-8

from Bio.Seq import Seq
from Bio import SeqIO
from fuzzysearch import find_near_matches
import statistics
import logging
import argparse
import os

# todo
# make output folder


# For each read in a fastq file
# find the positions of the primer, could be in forward or reverse orientation
# print the number of positions found and the lengths of the resulting reads, remember the inital one if not at start

parser = argparse.ArgumentParser(description='Please provide required inputs:')
required_group = parser.add_argument_group("Required Arguments")
required_group.add_argument("-f", "--fastq", required=True,  help="The full path to the fastq file (not zipped) to process. [REQUIRED]")
required_group.add_argument("-p", "--primer", required=True, help="The sequence of the forward primer used. [REQUIRED]")

optional_group = parser.add_argument_group("Optional Arguments")
optional_group.add_argument("-t", "--is_dumbell", action='store_true',
                            help="Is this a dumbell amplicon, i.e., will it forward and comp in same sequence.")
optional_group.add_argument("-c", "--no_dumbell_comp", action='store_true',
                            help="If this a dumbell amplicon, should the complemtary sequence NOT be complemented?")
optional_group.add_argument("-d", "--distance", default = 2, type = int,
                            help="Maximum Levenshtein distance for primer search [Default: 2]")
optional_group.add_argument("-mn", "--min_len", default = 2, type = int,
                            help="Minimum length of cut sequence [Default: 2]")
optional_group.add_argument("-mx", "--max_len", default = 9999999, type = int,
                            help="Maximum length of cut sequence [Default: 9999999]")
optional_group.add_argument("-o", "--out", help="Output folder [Default: cwd]")

args = parser.parse_args()
fastq_in = args.fastq
primer_seq = args.primer
primer_seq_c = Seq(primer_seq).reverse_complement() #also .complement()
is_dumbell = args.is_dumbell
dumbell_comp = args.no_dumbell_comp
max_d = args.distance
min_len = args.min_len
max_len = args.max_len
out = args.out

# fastq_in = "minIon_trial2.fastq" #"minIon.fastq"
# is_dumbell = False
# dumbell_comp = False
# primer_seq = "TCCTCCTCCGTTGTTGTTGTTG"
# primer_seq_c = Seq(primer_seq).complement() #also .reverse_complement()
# max_d = 2
# min_len = 20
# max_len = 5000

if out != None:
    if not os.path.isdir(out):
        try:
            os.makedirs(out, exist_ok=True)
            os.chdir(out)
        except:
            print("Could not make dir %s" % (out))

get_strt = lambda x: str(x).split(', ')[0].split('=')[1] # same as str(matches[1]).split(', ')[0].split('=')[1]
read_count = 0
primer_starts = []
passed_lengths = []
for record in SeqIO.parse(fastq_in, "fastq"):
    read_count += 1
    matches = find_near_matches(primer_seq, record.seq, max_l_dist=max_d)
    matches_comp = find_near_matches(primer_seq_c, record.seq, max_l_dist=max_d)
    if not is_dumbell:
        if len(matches) > 0 and len(matches_comp) > 0:
            if len(matches) >= len(matches_comp):
                logging.info("Found forward and comp primers for read %d, using highest which is forward" % (read_count))
            else:
                logging.info("Found forward and comp primers for read %d, using highest which is complement" % (read_count))
                matches = matches_comp
        elif len(matches) < 2:
            if len(matches_comp) >= 2: # means its comp strand sequenced
                matches = matches_comp
            else:
                if len(matches) == 1 or len(matches_comp) == 1: # means there is only one read
                    primer_starts.append(1)
                    SeqIO.write(record, "read_" + str(read_count) + "_single.fq", "fastq") # this does not append
                else:
                    primer_starts.append(0)
                pass
    else:
        db_order = [True] * len(matches) + [False] * len(matches_comp)
        matches = matches + matches_comp
        if len(matches) < 2:
            if len(matches) == 1 : # means there is only one read
                primer_starts.append(1)
                SeqIO.write(record, "read_" + str(read_count) + "_single.fq", "fastq") # this does not append
            else:
                primer_starts.append(0)
            pass
        else:
            if db_order[-1]: # add last one
                db_order.append(False)
            else:
                db_order.append(True)
    
    if len(matches) > 2:
        # So now matches has all the info
        strts = [0] + list(map(int, list(map(get_strt, matches)))) + [len(record.seq)]
        if not(dumbell_comp):
            [x for _,x in sorted(zip(sorted(strts[1:-1]), db_order))] # sort list based on other list

        strts.sort() # for the joining in dumbell method

        # lets get some stats
        #     rem duplicates if there are because of a perfect primer start and the addition of 0 above
        strts = list(dict.fromkeys(strts))
#         also need to do this for 
        if not(dumbell_comp):
            rem = [idx for idx, item in enumerate(strts[1:]) if item in strts[1:][:idx]]
#             db_order.pop(rem) #never used pop before, it prints the one it removes, easier than db_order[x:] + db_order[x+1:]
            db_order = [j for i, j in enumerate(db_order) if i not in rem]
        read_lngs = [j-i for i, j in zip(strts[:-1], strts[1:])]
        #maybe do mean, sd and median IRQ.
        logging.info("Mean lenth for read %d is %d, with a min of %d, a max of %d before filtering (if user selected)" % (read_count, statistics.mean(read_lngs), min(read_lngs), max(read_lngs)))

        # filter
        if sum(1 for x in read_lngs if min_len <= x <= max_len) < 2:
            primer_starts.append(1)
            ext = "_single.fq"
        else:
            ext = "_multiple.fq"

        records_tmp = []
        passed_lengths_ind = []
        for i in range(0, len(strts)-1):
            cut_len = strts[i+1] - strts[i]
            if min_len <= cut_len <= max_len:
                if not(dumbell_comp) and db_order[i]: # we need to the complement those reads with the complement of the primer sequence (qual score will stay the same)
    #                 Figure out which ones they are and do this
                    record.seq = record.seq.reverse_complement()
                records_tmp.append(record[strts[i]:strts[i+1]])
                passed_lengths_ind.append(cut_len)

        passed_lengths.append(passed_lengths_ind)        
        SeqIO.write(records_tmp, "read_" + str(read_count) + ext, "fastq") # this does not append
        # this can now be passed to the other script to make a consensus
        primer_starts.append(len(strts)-1) #or -2?

try:
    diff_lens = [max(i)-min(i) for i in passed_lengths]
    logging.info("The median distance between the shortest and longest read lengths in the cut sequences is %d" % statistics.median(diff_lens))
    with_seq = len([i for i in primer_starts if i > 1]) 
    print("From the %d input reads, %d had more than one primer start and there was a media of %d primer starts. The mean length of the median cut reads per sample output is %d" % (read_count, with_seq, statistics.mean(primer_starts), statistics.median([statistics.median(i) for i in passed_lengths])))
except:
    print("cant calc median, probaly no data")


