#!/bin/python
import gzip
import argparse
import re
import os
#import numpy as np

#Set input variables
def get_args():
    parser = argparse.ArgumentParser(description="Identify reads that map to the same chromosome and position of a genome and remove duplicated contigs")
    parser.add_argument("-f", "--file", help="required arg, absolute file path", type=str)
    parser.add_argument("-u", "--umi", help="optional arg, designates file containing the list of UMIs (unset if randomers instead of UMIs)", type=str)
#    parser.add_argument("-o", "--out_file", help="Desired outfile path + name", type=str)
    parser.add_argument("-p", "--paired", help="optional arg, designates file is paired end (not single-end) / Input options: paired, single", type=str)
    parser.add_argument("-t", "--sort_val", help="Is the sam file sorted? (1 = yes | 0 = No) *Make sure samtools module is loaded if sam file is not sorted*", type=int)
#    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
    return parser.parse_args()

args = get_args()
# Setting user inputted variables
Sam_file = args.file
index_file = args.umi
Sorted_Val = args.sort_val
Paired = args.paired
#out_file = args.out_file

# Format outfile name
file_name = Sam_file.split(".sam")[0]
out_file = file_name + "_deduped.sam"

# Check if paired functionality is requested
if Paired == "paired" or Paired == "Paired":
    print("No paired-end functionality, sorry!")
    exit()

# Sort inputted sam file if desired
if Sorted_Val == 0:
    # Samtool commands:
    Sam_2_Bam = "samtools view -S -b " + Sam_file + " > temp.bam" # convert sam to bam file
    Sort_Bam = "samtools sort temp.bam -O sam -o sorted.sam"
    os.system(Sam_2_Bam)
    os.system(Sort_Bam)

#A set to hold the known UMIs
Known_UMIs = set()

# Store all the known UMIs into a set
with open(index_file) as umi:
    for line in umi:
        line = line.strip()
        Known_UMIs.add(line)

def get_tuple(record):
    '''Take in a record and create a tuple containing the UMI and 5' Position'''
    # Define directionality variable
    Dir = 0
    # Tuple to store UMI and Adjusted Left Most Position
    tup = tuple()

    # Get UMI 0
    headerline = record[0].split(":")
    UMI = headerline[7]

    # Check if UMI is valid
    if UMI in Known_UMIs:

        # Get the directionality
        bit_flag = int(record[1])
        if ((bit_flag & 16) == 16):
            Dir = 1

        # Get left most aligned position
        left_pos = int(record[3])

        # Get the corrected 5' starting position
        cigar = record[5]
        cig_split = re.findall(r'\d+[MIDNSHP=X]', cigar)

        # Adjust for soft clipping on left side if direction is forward (1)
        if Dir == 0:
            if 'S' in cig_split[0]:
                left_pos = left_pos - int(cig_split[0][:-1])

        # Adjust for soft clipping, insertions, and deletions if direction if reverse
        elif Dir == 1:
            if 'S' in cig_split[0]:
                cig_split.pop(0)
            for obj in cig_split:
                if 'I' not in obj:
                    left_pos += int(obj[:-1])

        # Create tuple to be returned
        tup = (UMI, left_pos)

    # If UMI is not in the known UMI set, make Dir = -1
    else:
        Dir = -1

    #Return the tuple and the directionality
    return (tup, Dir)

# Setting Global Variables
# Keeps track of the previous chromosome value in order to flush the record storing sets
prev_chromosome = 0
# Used to identify the directionality of the reads.
Forward = 0
Reverse = 1
# Sets used to keep track of records that have been written out to the out file
Forward_set = set()
Reverse_set = set()

# Open file to write to
write_file = open(out_file, "a")

with open("sorted.sam") as body:
    for line in body:
        split_line = line.split()

        # Write headerlines to out output file
        if '@' in split_line[0]:
            write_file.write(line)

        else:
            # Parse out the chromosome value and compare it to stored chromosome value
            new_chromosome = split_line[2]

            # Compare if chromosome values are the same or different. Sets are flushed to reduce memory storage.
            if new_chromosome != prev_chromosome:
                prev_chromosome = new_chromosome
                Forward_set.clear()
                Reverse_set.clear()

            # Obtain the tuple used to check if recorord exists and the directionality of the aligned read.
            record_tup, direction = get_tuple(split_line)

            # Do not write records to file if UMI is not valid
            if direction != -1:
                # Write record to the output file if it is not a PCR duplicate
                if direction == Forward and record_tup not in Forward_set:
                    write_file.write(line)
                    Forward_set.add(record_tup)
                if direction == Reverse and record_tup not in Reverse_set:
                    write_file.write(line)
                    Reverse_set.add(record_tup)

write_file.close()

# remove temporary files from current directory
os.remove("temp.bam")
os.remove("sorted.sam")
