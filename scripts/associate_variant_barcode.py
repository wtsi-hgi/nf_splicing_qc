#-- import modules --#
import argparse
import gzip
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from collections import Counter

#-- functions --#
def hamming_distance(str1, str2):
    """
    Calculate the Hamming distance between two strings.
    """
    if len(str1) != len(str2):
        return max(len(str1), len(str2))
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def match_approximate(seq, pattern, max_mismatches):
    """
    Find approximate match of pattern in seq allowing max_mismatches.
    Returns start index of match or -1 if not found.
    """
    k = len(pattern)
    for i in range(len(seq) - k + 1):
        window = seq[i:i+k]
        if hamming_distance(window, pattern) <= max_mismatches:
            return i
    return -1

def extract_sequence(seq, up_seq, down_seq, max_mismatches):
    """
    Extract substring from seq between approximate matches of up_seq and down_seq.
    """
    start_idx = match_approximate(seq, up_seq, max_mismatches)
    if start_idx == -1:
        return "upstream not found"
    start_idx += len(up_seq)

    end_idx = match_approximate(seq[start_idx:], down_seq, max_mismatches)
    if end_idx == -1:
        return "downstream not found"
    end_idx += start_idx

    return seq[start_idx:end_idx]

def reverse_complement(seq):
    """
    Generate the reverse complement of a DNA sequence.
    """
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def check_barcode(barcode_seq, barcode_temp, max_mismatches):
    """
    Check if barcode_seq matches the barcode_temp allowing max_mismatches
    """
    if len(barcode_seq) != len(barcode_temp):
        return False

    mismatch_count = 0
    for s_char, p_char in zip(barcode_seq, barcode_temp):
        if p_char != 'N':
            if s_char != p_char:
                mismatch_count += 1
                if mismatch_count > max_mismatches:
                    return False

    return True

def process_read_pair(read1, read2, variant_up, variant_down, barcode_up, barcode_down, max_mismatches):
    """
    Extract variant from read1 between variant_up and variant_down.
    Extract barcode from read2 between barcode_up and barcode_down.
    """
    seq1 = str(read1.seq)
    seq2 = str(read2.seq)

    variant_seq = extract_sequence(seq1, variant_up, variant_down, max_mismatches)
    barcode_seq = extract_sequence(seq2, barcode_up, barcode_down, max_mismatches)
    
    return (variant_seq, barcode_seq)

def read_fastq_in_chunk(r1_path, r2_path, variant_up, variant_down, barcode_up, barcode_down, max_mismatches, chunk_size, threads):
    """
    Read paired-end FASTQ files in chunks and process each pair of reads in parallel.
    """
    open_func = gzip.open if r1_path.endswith(".gz") else open

    with open_func(r1_path, "rt") as r1_handle, open_func(r2_path, "rt") as r2_handle:
        r1_iter = SeqIO.parse(r1_handle, "fastq")
        r2_iter = SeqIO.parse(r2_handle, "fastq")

        while True:
            r1_chunk, r2_chunk = [], []
            try:
                for _ in range(chunk_size):
                    r1_chunk.append(next(r1_iter))
                    r2_chunk.append(next(r2_iter))
            except StopIteration:
                pass

            if not r1_chunk or not r2_chunk:
                break

            with ThreadPoolExecutor(max_workers=threads) as executor:
                results = list(executor.map(
                    lambda args: process_read_pair(*args),
                    zip(r1_chunk, r2_chunk,
                        [variant_up]*len(r1_chunk),
                        [variant_down]*len(r1_chunk),
                        [barcode_up]*len(r1_chunk),
                        [barcode_down]*len(r1_chunk),
                        [max_mismatches]*len(r1_chunk))
                ))

            yield results

            if len(r1_chunk) < chunk_size:
                break

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Extract variant and barcode from paired-end FASTQ files.")
    parser.add_argument("--read1",            type = str, required = True, help = "Read 1 FASTQ file")
    parser.add_argument("--read2",            type = str, required = True, help = "Read 2 FASTQ file")
    parser.add_argument("--variant_len",      type = int, required = True, help = "Length of the variant sequence")
    parser.add_argument("--variant_up",       type = str, required = True, help = "Upstream flank for variant sequence")
    parser.add_argument("--variant_down",     type = str, required = True, help = "Downstream flank for variant sequence")
    parser.add_argument("--barcode_temp",     type = str, required = True, help = "Template for barcode sequence (e.g., 'NNATNNNNATNNNNATNNNN')")
    parser.add_argument("--barcode_up",       type = str, required = True, help = "Sequence before barcode in read2")
    parser.add_argument("--barcode_down",     type = str, required = True, help = "Sequence after barcode in read2")
    parser.add_argument("--barcode_check",    action="store_true",         help = "Enable barcode checking against template")
    parser.add_argument("--barcode_mismatch", type = int, default = 1,     help = "Number of mismatches allowed in barcode checking")
    parser.add_argument("--max_mismatches",   type = int, default = 2,     help = "Max mismatches allowed in up/down matches")
    parser.add_argument("--min_barcov",       type = int, default = 1,     help = "Minimum coverage for barcode-variant association")
    parser.add_argument("--chunk_size",       type = int, default = 5000,  help = "Chunk size for processing reads")
    parser.add_argument("--threads",          type = int, default = 12,    help = "Number of threads")
    parser.add_argument("--output_file",      type = str, default = "variant_barcode_association.txt", help = "Output file")
    args = parser.parse_args()

    all_results = []
    for i, result_chunk in enumerate(read_fastq_in_chunk(args.read1, args.read2, 
                                                         args.variant_up.upper(), args.variant_down.upper(), 
                                                         args.barcode_up.upper(), args.barcode_down.upper(), 
                                                         args.max_mismatches, args.chunk_size, args.threads)):
        print(f"Processed chunk {i+1} with {len(result_chunk)} read pairs")
        all_results.extend(result_chunk)

    count_varup_notfound = 0
    count_vardown_notfound = 0
    count_barup_notfound = 0
    count_bardown_notfound = 0
    count_variant_length = 0
    if args.barcode_check:
        count_barcode_pattern = 0
    count_effective_reads = 0
    all_results_filtered = []
    for (variant_seq, barcode_seq) in all_results:
        if variant_seq == "upstream not found":
            count_varup_notfound += 1
        elif variant_seq == "downstream not found":
            count_vardown_notfound += 1
        else:
            if barcode_seq == "upstream not found":
                count_barup_notfound += 1
            elif barcode_seq == "downstream not found":
                count_bardown_notfound += 1
            else:
                if len(variant_seq) != args.variant_len:
                    count_variant_length += 1
                else:
                    if args.barcode_check:
                        if check_barcode(reverse_complement(barcode_seq), args.barcode_temp.upper(), args.barcode_mismatch):
                            count_effective_reads += 1
                            all_results_filtered.append((variant_seq, reverse_complement(barcode_seq)))
                        else:
                            count_barcode_pattern += 1
                    else:
                        count_effective_reads += 1
                        all_results_filtered.append((variant_seq, reverse_complement(barcode_seq)))

    print(f"Total reads processed: {len(all_results)}")
    print(f"Total reads with variant upstream not found: {count_varup_notfound}")
    print(f"Total reads with variant downstream not found: {count_vardown_notfound}")
    print(f"Total reads with barcode upstream not found: {count_barup_notfound}")
    print(f"Total reads with barcode downstream not found: {count_bardown_notfound}")
    print(f"Total reads with variant length inconsistent: {count_variant_length}")
    if args.barcode_check:
        print(f"Total reads with barcode pattern mismatch: {count_barcode_pattern}")
    print(f"Total effective reads: {count_effective_reads}")

    variant_barcode_association = Counter(all_results_filtered)

    output_file = args.output_file
    with open(output_file, "w") as f:
        f.write("Barcode\tVariant\tCount\n")
        for (variant, barcode), count in variant_barcode_association.items():
            f.write(f"{barcode}\t{variant}\t{count}\n")
