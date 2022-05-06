#!/usr/bin/env python3

import argparse
import Bioinfo
import re

################################## ARGPARSE ##################################
def get_args():
    parser = argparse.ArgumentParser(description='Remove PCR duplicates from SAM file and output new SAM file.')
    parser.add_argument('-f', '--file', type=str, required=True, 
                        help='absolute path to SAM file')
    parser.add_argument('-u', '--umi', type=str, required=True, 
                        help='file containing the list of UMIs (unset if randomers instead of UMIs)')
    parser.add_argument('-p', '--paired', required=False,  action="store_true",
                        help='flag if file is paired end')
    return parser.parse_args()

# Assign variables from arguments in the command line
args = get_args()
sam_file = args.file
umi_file = args.umi
output_file = sam_file.split('.sam')[0].split('/')[-1] + "_deduped.sam"
if args.paired == True:
    print("Error: This script may only be run with single-end data.")
    quit()

################################## UMIs ##################################

## open umi file, store umis as set, close file. 
f = open(umi_file, "r")
umi_set = set()
for line in f:
    line = line.rstrip()
    umi_set.add(line)
f.close()

################################## COMPARE SAMLINES ##################################

## dictionary for comparing reads at the same position
## check_dup_dict: key="chr POS revc UMI", value=sam line
no_dup_dict = {}
chr_tracker = ""  # keeps track of genome position as we go through chromosome

## open input file
fh_input = open(sam_file, "r")
fh_output = open (output_file, "a")
fh_output.truncate(0) # wipe output file if it already exists

## diagnostic numbers
num_headers = 0     # number of header lines
umi_err = 0         # number of reads with bad umis
unique_reads = 0    # number of unique reads
dup_reads = 0       # number of duplicated reads

for line in fh_input:
    if line.startswith("@"): # if header line, write directly to output file
        fh_output.write(line)
        num_headers += 1
    else:
        line_dict = Bioinfo.dice_head(line)
        if line_dict['umi'] not in umi_set: # skip line if umi is not in known umi set
            umi_err += 1
        else:
            if line_dict['chr'] == chr_tracker: # if this read is within the same chromosome as previous
                # convert POS to real position and store in dictionary
                new_pos = Bioinfo.to_true_pos(line_dict['POS'], line_dict['CIGAR'], line_dict['revc'])
                line_dict['POS'] = str(new_pos) # replaces position in line_dict with true read position
                no_dup_key = ' '.join([str(i) for i in [line_dict['chr'], line_dict['POS'], line_dict['revc'], line_dict['umi']]])
                if no_dup_key in no_dup_dict:
                    dup_reads += 1
                else:
                    # make a new entry if it does not exist
                    no_dup_dict[no_dup_key] = line
                    # write to file the first unique line
                    fh_output.write(line)
                    unique_reads += 1
            else: # this read is in a different chromosome
                chr_tracker = line_dict['chr'] # update pos_tracker value to match new position
                # reset dictionary and add new line for first entry
                no_dup_dict = {}
                # convert POS to real position and store in dictionary
                new_pos = Bioinfo.to_true_pos(line_dict['POS'], line_dict['CIGAR'], line_dict['revc'])
                line_dict['POS'] = str(new_pos) # replaces position in line_dict with true read position
                no_dup_key = ' '.join([str(i) for i in [line_dict['chr'], line_dict['POS'], line_dict['revc'], line_dict['umi']]])
                # write to file the first unique line
                fh_output.write(line)
                unique_reads += 1
fh_input.close()
fh_output.close()

## print run information
print("Input file:", sam_file.split('/')[-1])
print("Output file:", output_file)
print("Number of header lines:", num_headers)
print("Number of unique reads (final output):", unique_reads)
print("Number of reads with UMI errors:", umi_err)
print("Number of duplicate reads removed:", dup_reads)

