#!/usr/bin/env python3

import numpy as np
import re
from itertools import combinations

# Constants:
DNA_bases = ("A", "T", "G", "C")
RNA_bases = ("A", "U", "G", "C")
test_seq = "TGCAGGTTGAGTTCTCGCTGTCCCGCCCCCATCTCTTCTCTTCCAGTCTGGCTCTGGAGCAGTTGAGCCCAGCTCAGGTCCGCCGGAGGAGACCG"
test_phred_score = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"
short_seq = "TGCAGGTT"
bar_dict = {"B1": "GTAGCGTA", "A5": "CGATCGAT", "C1": "GATCAAGG", "B9": "AACAGCGA"}
sam_line = "NS500451:154:HWKTMBGXX:1:11101:24260:1121:GATCAAGG	0	2	76814284	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU"
umi_set = {v for v in bar_dict.values()}

#####################################################################################

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter)-33

if __name__ == "__main__":
	assert convert_phred("A") == 32, "Phred score incorrect"
	assert convert_phred("@") == 31, "Phred score incorrect"
	assert convert_phred("#") == 2, "Phred score incorrect"
	print("PASS\tconvert_phred")

#####################################################################################

def qual_score(phred_score: str) -> float:
    '''Converts the phred score line from a read and calculates average quality score of read'''
    psum = 0
    plen = len(phred_score)
    for i in range(plen):
        psum += convert_phred(phred_score[i])
    av_q = psum/plen
    return av_q

if __name__ == "__main__":
	assert int(qual_score(test_phred_score)) == 37, "Average qual score incorrect"
	print("PASS\tqual_score")

#####################################################################################

def validate_base_seq(seq: str, RNAflag: bool = False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = set(seq.upper())
    model = set("AUGC" if RNAflag else "ATGC")
    return seq <= model  # if seq is within model, will be True

if __name__ == "__main__":
	assert validate_base_seq(
	    "AATAGAT") == True, "Validate base seq does not work on DNA"
	assert validate_base_seq(
	    "AAUAGAU", True) == True, "Validate base seq does not work on RNA"
	print("PASS\tvalidate_base_seq")

#####################################################################################

def gc_content(seq: str) -> float:
    '''Calculates the GC content of a DNA sequence'''
    seq = seq.upper()
    gc = seq.count("C") + seq.count("G")
    gc_cont = gc/len(seq)
    return gc_cont

if __name__ == "__main__":
	assert gc_content("GCGCGC") == 1
	assert gc_content("AATTATA") == 0
	assert gc_content("GCATGCAT") == 0.5
	print("PASS\tgc_content")

#####################################################################################

def oneline_fasta(infile: str, outfile: str):
    '''
    Converts input fasta file into nicer format fasta.
    Works by only adding newlines right before each header
    unless header is very first line in input file.
    '''
    ofh = open(outfile, "a")
    ofh.truncate(0)  # clears file before appending
    with open(infile, "r") as ifh:
        i = 0
        for line in ifh:
            i += 1
            line = line.strip("\n")
            if line.startswith(">"):
                if i != 1:
                    ofh.write("\n")
                ofh.write(line+"\n")
            else:
                ofh.write(line)
    ofh.close()
    print(infile, "converted to", outfile)

#####################################################################################
def read_stats(file: str) -> tuple:
    '''Returns average read length of input fastq file.'''
    with open(file, "r") as f:
        nt_count, max_rlen, min_rlen = 0, 0, 0
        i = 0
        for line in f:
            i += 1
            if i==2:
                min_rlen = len(line.strip()) # set a number for min len
            if i % 4 == 2: # if it's a read seq line
                read_len = len(line.strip())
                nt_count += read_len # add nucleotide count to counter
                if read_len > max_rlen:
                    max_rlen = read_len
                if read_len < min_rlen:
                    min_rlen = read_len
        readnum = i/4 # number of reads in file
        av_rlen = nt_count/readnum
        rstat_tup = (av_rlen, min_rlen, max_rlen)
        return rstat_tup

if __name__ == "__main__":
	assert read_stats("test.fastq") == (319, 214, 402)
	print("PASS\tread_stats")

#####################################################################################

def revc (seq: str) -> str:
	'''converts DNA seq to reverse complement'''
	seq=seq.upper()
	nucd = {"A":  "T", "T": "A", "C": "G", "G": "C", "N": "N"}
	rseq=seq.translate(str.maketrans(nucd))[::-1] 
	return rseq

if __name__ == "__main__":
	assert revc(short_seq) == "AACCTGCA"
	print("PASS\trevc")

#####################################################################################

def ham_dist(seq1: str, seq2: str) -> int:
    '''calculates hamming dist between two sequences'''
    dict_hdist = sum(s1!=s2 for s1, s2 in zip(seq1, seq2))
    return dict_hdist

if __name__ == "__main__":
	assert ham_dist(short_seq, "TGCAAAAA") == 4, "Hamming distance incorrect"
	print("PASS\tham_dist")

#####################################################################################

def min_hdist(seq_list: list) -> int:
    '''calculates minimum hamming dist between sequences in a list'''
    combo_list = combinations(seq_list, 2)
    minh = len(seq_list[1]) # initialize min ham dist w/length of an index
    for pair in combo_list:
        h = ham_dist(pair[0], pair[1])
        if h < minh:
            minh = h
    return minh

if __name__ == "__main__":
    bars = [b for b in bar_dict.values()]
    assert min_hdist(bars) == 5, "Hamming distance incorrect"
    print("PASS\tmin_hdist")

#####################################################################################

def bc_correct(qbar_seq: str, barcode_lst: list) -> str:
    '''corrects sequence of barcode given a min hamming dist and barcode dict'''
    hdist = min_hdist(barcode_lst) # min ham dist in barcode list
    corrected_index = qbar_seq # by default returns input barcode if no correction
    # calculate & store min hamming dist between qbar_seq and each barcode in list
    # hdict: key=barcode, value=hdist
    hdict = {}
    for b in barcode_lst:
        h = ham_dist(qbar_seq, b)
        hdict[b] = h
    minham = min(x for x in hdict.values())
    # if this ham dist is < user-provided min ham dist, find best match
    if minham < hdist:
        # make sure there is only one best barcode, else do no correction
        minham_count = min(v for v in hdict.values())
        if minham_count == 1:
            corrected_index = [k for k,v in hdict.items() if v == minham][0]
    return corrected_index

if __name__ == "__main__":
    bars = [b for b in bar_dict.values()]
    assert bc_correct("GTAGCGTT", bars) == "GTAGCGTA", "Phred score incorrect"
    print("PASS\tbc_correct")

#####################################################################################

def dice_head (sam_head: str) -> dict:
    '''converts SAM line to list of components needed for analysis'''
    sam_head.count('\t')
    [QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL,misc] = sam_head.split("\t", maxsplit=11)
    UMI = QNAME.split(':')[7]
    FLAG, POS = [ int(x) for x in [FLAG, POS]]
    if ((FLAG & 4) == 4):
        unmapped = True
    else:
        unmapped = False
    if ((FLAG & 16) == 16): 
        revc = True
    else: 
        revc = False
    return {"umi":UMI, "FLAG":FLAG, "chr":RNAME, "POS":POS, "CIGAR":CIGAR, "unmapped":unmapped, "revc":revc}

if __name__ == "__main__":
    assert dice_head(sam_line) == {'umi': 'GATCAAGG', 'FLAG': 0, 'chr': '2', 'POS': 76814284, 'CIGAR': '71M', 'unmapped': False, 'revc': False}, "Incorrect values in list"
    sam_line = "NS500451:154:HWKTMBGXX:1:11101:24260:1121:GATCAATG	18	2	76814284	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU"
    assert dice_head(sam_line) == {'umi': 'GATCAATG', 'FLAG': 18, 'chr': '2', 'POS': 76814284, 'CIGAR': '71M', 'unmapped': False, 'revc': True}, "Incorrect values in list"
    print("PASS\tdice_head")

#####################################################################################

def to_right_pos (sam_pos: int, cigar: str) -> int:
    '''Converts SAM position with soft clipping (from CIGAR string) to right (genomic 3') position'''
    cigar_letters = re.findall(r"\D", cigar)
    cigar_values = re.findall(r"\d+", cigar)
    total_M = 0
    total_I = 0
    total_N = 0
    total_D = 0
    total_S = 0
    left_S = 0
    for i in range(len(cigar_letters)):
        if cigar_letters[i] == "M" or cigar_letters[i] == "X" or cigar_letters[i] == "=":
            total_M = total_M + int(cigar_values[i])
        elif cigar_letters[i] == "I":
            total_I = total_I + int(cigar_values[i])
        elif cigar_letters[i] == "N":
            total_N = total_N + int(cigar_values[i])
        elif cigar_letters[i] == "D":
            total_D = total_D + int(cigar_values[i])
        elif cigar_letters[i] == "S":
            total_S = total_S + int(cigar_values[i])
            if i == 0:
                left_S = int(cigar_values[0])
    read_length = total_M + total_I + total_S
    right_pos = sam_pos + read_length - 1 + total_D - left_S + total_N - total_I
    #right_pos = sam_pos + total_M + total_S - left_S + total_D + total_N - 1
    return right_pos

#left_pos = 2
#cigar = "3S10M5S"
#print(to_right_pos (left_pos, cigar))

if __name__ == "__main__":
    left_pos = 2
    cigar = "8M"
    assert to_right_pos (left_pos, cigar) == 9, "Incorrect position value (M)"
    cigar = "1S6M"
    assert to_right_pos (left_pos, cigar) == 7, "Incorrect position value (left S)"
    cigar = "6M3S"
    assert to_right_pos (left_pos, cigar) == 10, "Incorrect position value (right S)"
    cigar = "50M10N5M"
    assert to_right_pos (left_pos, cigar) == 66, "Incorrect position value (N)"
    cigar = "3M2I3M"
    assert to_right_pos (left_pos, cigar) == 7, "Incorrect position value (I)"
    cigar = "3M2D3M"
    assert to_right_pos (left_pos, cigar) == 9, "Incorrect position value (D)"
    cigar = "3S10M5S"
    assert to_right_pos (left_pos, cigar) == 16, "Incorrect position value (left & right S)"
    print("PASS\tto_right_pos")

#####################################################################################

def to_true_pos (soft_pos: int, cigar: str, revc: bool) -> int:
    '''Converts position with soft clipping (from CIGAR string) to true position'''
    cigar_letters = re.findall(r"\D", cigar)
    cigar_values = re.findall(r"\d+", cigar)
    if revc == False: # top strand
        if cigar_letters[0] == 'S':
            true_pos = soft_pos - int(cigar_values[0])
            return true_pos
        else:
            true_pos = soft_pos
            return true_pos
    else: # bottom strand
        true_pos = to_right_pos (soft_pos, cigar)
        return true_pos

#sam_pos = 50
#cigar = "3S15M"
#print(to_true_pos(sam_pos, cigar, True))

if __name__ == "__main__":
    sam_pos = 50
    cigar = "18M"
    assert to_true_pos(sam_pos, cigar, False) == 50, "Incorrect position value"
    cigar = "2S18M"
    assert to_true_pos(sam_pos, cigar, False) == 48, "Incorrect position value"
    cigar = "3S15M"
    assert to_true_pos(sam_pos, cigar, True) == 64, "Incorrect position value (left S), bottom strand"
    cigar = "3S10M5S"
    assert to_true_pos(sam_pos, cigar, True) == 64, "Incorrect position value (two S), bottom strand"
    print("PASS\tto_true_pos")

#####################################################################################

