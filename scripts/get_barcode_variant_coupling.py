#import the modules
import sys
import pandas as pd
import numpy as np
import io
import os
from collections import Counter
import argparse
import pickle
import gzip

# https://github.com/regev-lab/interpretable-splicing-model/blob/main/data_preprocessing/utils.py

#define a function that measure the hamming distance between 2 strings of characters
def hamming(str1, str2):
    if len(str1) != len(str2):
        return max(len(str1), len(str2))
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

#define a function that make the rev complement (necessary to get the bc in the + strand orientation)
def revcomp(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

#define variables
filenames = [sys.argv[1], sys.argv[2]] #file to be analysed (already specified in the .sh file)
barcodeexondict = {} #empty variable to store good barcode
badbarcodes = [] #empty variable to store bad barcode
counter = 0 # initialize a variable to 0 it will be used to count how many times a specific bc-var association is seen
maxerrs = 2 #maximum number of mismatch allowd


#Variant (read1)
upstreamflank = sys.argv[3] #constant sequence before the start of the variable sequence (already specified in the .sh file)
downstreamflank = sys.argv[4] #constant sequence after the end of the variable sequence (already specified in the .sh file)
exonposition = len(upstreamflank) #where is the start of the variant sequence (1 nt after the end of the constant part rigth upstream)
# exonlen = sys.argv[5] #length of the variant (already specified in the .sh file)
# lenfile1 = exonposition+exonlen+len(downstreamflank)


#Barcode (read2)
barcodeflank="AAGCTTGCATCGAATCAGTAG" #constant sequence after the end of the barcode' 
exon7="TTCCTTTCTGTGCTTTCTGC" #constant sequence before the start of the barcode' 
barcodelen = 38 #length of the barcode
lenfile2 = len(exon7)+barcodelen+len(barcodeflank) # length of the sequence with the barcode + 2 constant parts (barcodeflank, exon7)



#initialize a series of variable to keep track of how many good reads we are keeping + how many bad reasd we're discrading ad why (save in the output file)
readsused = 0
readstossed = 0
# Counters for different reasons for tossing reads
tossed_due_to_length_file2 = 0 #discared because bc has an incorrect length 
tossed_due_to_barcodeflank_mismatch = 0 #discared for too many mismatch in constant sequence after the end of the barcode
tossed_due_to_exon7_mismatch = 0 #discared for too many mismatch in constant sequence before the start of the barcode
tossed_due_to_length_file1 = 0 #discared for incorrect length of the variant
tossed_due_to_upstreamflank_mismatch = 0 #discared for too many mismatch in constant sequence before the start of the variant
tossed_due_to_downstreamflank_mismatch = 0 #discared for too many mismatch in constant sequence after the end of the variant

#START of the actual code
# open the files
def gen_line(filename):
    with gzip.open(filename, 'rb') as f:
        for line in f:
            line = line.decode("utf-8")
            yield line.strip()


#number of line in the file
gens = [gen_line(n) for n in filenames]

#process each line independently (one after the other)
for file1_line, file2_line in zip(*gens):
    if counter % 10e6 == 0: #every 1M line processed update the variable define above to count the number of good/bad reads
        print(f"processed {counter} lines")
        print(f"Total reads used: {readsused}")
        print(f"Total reads tossed: {readstossed}")
        print(f"Reasons for tossing reads:")
        print(f" - Incorrect length in file2: {tossed_due_to_length_file2}")
        print(f" - Barcodeflank mismatch: {tossed_due_to_barcodeflank_mismatch}")
        print(f" - Exon7 mismatch: {tossed_due_to_exon7_mismatch}")
        print(f" - Incorrect length in file1: {tossed_due_to_length_file1}")
        print(f" - Upstreamflank mismatch: {tossed_due_to_upstreamflank_mismatch}")
        print(f" - Downstreamflank mismatch: {tossed_due_to_downstreamflank_mismatch}")
        with open("lib9_barcode_dict_cluster.pickle","wb") as f:
            pickle.dump(barcodeexondict, 
                        f) 
    counter += 1 #

    exon = None
    barcode = None
    length = None
    
    # first get the barcode
    #check the bc has the correct length
    if (len(file2_line[(len(barcodeflank)):(len(barcodeflank) + barcodelen)]) !=barcodelen): 
        tossed_due_to_length_file2 += 1
        continue

    #check constant part has not +2 mismatch
    if hamming(file2_line[:len(barcodeflank)], barcodeflank) > maxerrs: 
        readstossed += 1
        tossed_due_to_barcodeflank_mismatch += 1
        continue
    
    #check constant part has not +2 mismatch
    if hamming(file2_line[(len(barcodeflank) + barcodelen):(len(barcodeflank) + barcodelen + len(exon7))], exon7) > maxerrs: 
        readstossed += 1
        tossed_due_to_exon7_mismatch += 1
        continue
    
    #check the var has the correct length
    # if (len(file1_line[exonposition :(exonposition + exonlen)]) < lenfile1-5 and len(file1_line[exonposition :(exonposition + exonlen)]) > lenfile1+5): 
    #     tossed_due_to_length_file1 += 1
    #     continue

    #check constant part has not +2 mismatch
    if hamming(file1_line[:len(upstreamflank)], upstreamflank) > maxerrs: 
        tossed_due_to_upstreamflank_mismatch += 1
        readstossed += 1
        continue

    #check constant part has not +2 mismatch
    if hamming(file1_line[(len(file1_line)-len(downstreamflank)):], downstreamflank) > maxerrs: 
        tossed_due_to_downstreamflank_mismatch += 1
        readstossed += 1
        continue
    
    exon = file1_line[exonposition :(len(file1_line)-len(downstreamflank))] #store the variant seq
    rcbarcode = file2_line[(len(barcodeflank)):(len(barcodeflank) + barcodelen)] #store the bc seq
    barcode = revcomp(rcbarcode) #store the rev_comp of the bc seq

    if not barcode in barcodeexondict.keys():
        barcodeexondict[barcode] = [
            Counter(),
            0,
        ]  # second coordinate counts number of bad reads

    barcodeexondict[barcode][0][exon] += 1
    readsused += 1

# Final summary
print(f"Total reads processed: {counter}")
print(f"Total reads used: {readsused}")
print(f"Total reads tossed: {readstossed}")
print(f"Reasons for tossing reads:")
print(f" - Incorrect length in file2: {tossed_due_to_length_file2}")
print(f" - Barcodeflank mismatch: {tossed_due_to_barcodeflank_mismatch}")
print(f" - Exon7 mismatch: {tossed_due_to_exon7_mismatch}")
print(f" - Incorrect length in file1: {tossed_due_to_length_file1}")
print(f" - Upstreamflank mismatch: {tossed_due_to_upstreamflank_mismatch}")
print(f" - Downstreamflank mismatch: {tossed_due_to_downstreamflank_mismatch}")



barcodes=list(barcodeexondict.keys())

print("total number of barcodes:", len(barcodes))

barthresh= sys.argv[6] #num min reads to keep a bc-var association

barcodeskeep=[]
exonskeep=[]
counts=[]

#this keep the most common association per each bc-var
for i in range(len(barcodes)):
    curcount=barcodeexondict[barcodes[i]][0].most_common(1)[0][1]
    if curcount >= barthresh:
        barcodeskeep.append(barcodes[i])
        exonskeep.append(barcodeexondict[barcodes[i]][0].most_common(1)[0][0])
        counts.append(curcount)

print("total number of barcodes with at least 10 read:", len(barcodeskeep))


bedf=pd.DataFrame({"barcode" : barcodeskeep,
                   "exon" : exonskeep,
                   "count" : counts})

#save the table with barcode seq - variant seq - # reads
bedf.to_csv(sys.argv[7]+"_barcode_exon_pairs_min"+sys.argv[7]+"reads_.txt",index=False,sep="\t")


