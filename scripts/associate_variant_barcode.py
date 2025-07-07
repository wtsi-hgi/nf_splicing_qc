#-- import modules --#
import io
import sys
import argparse
import subprocess
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from collections import Counter

#-- functions --#
def pigz_open(path: str):
    """
    Open a gzip file using pigz for faster decompression
    Parameters:
        -- path: file path to the gzip file
    Returns:
        -- io_wrapper: a TextIOWrapper for reading the decompressed file
    """
    return subprocess.Popen(["pigz", "-dc", path], stdout = subprocess.PIPE)

def fastq_iter(handle):
    """
    FASTQ parser yielding (header, seq, qual)
    Parameters:
        -- handle: file handle for the FASTQ file
    Yields:
        -- (header, seq, qual): a tuple containing the header, sequence, and quality
    """
    while True:
        try:
            header = next(handle).strip()
            seq    = next(handle).strip()
            next(handle)  # skip '+'
            qual   = next(handle).strip()
            yield (header, seq, qual)
        except StopIteration:
            break

def hamming_distance(str1: str, str2: str) -> int:
    """
    Calculate the Hamming distance between two strings
    Parameters:
        -- str1: first string
        -- str2: second string
    Returns:
        -- int: the Hamming distance, or the maximum length if they differ in length
    """
    if len(str1) != len(str2):
        return max(len(str1), len(str2))
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def match_approximate(seq: str, pattern: str, max_mismatches: int) -> int:
    """
    Find approximate match of pattern in seq allowing max_mismatches.
    Returns start index of match or -1 if not found.
    Parameters:
        -- seq: the sequence to search in
        -- pattern: the pattern to match
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- int: start index of the match or -1 if not found
    """
    k = len(pattern)
    for i in range(len(seq) - k + 1):
        window = seq[i:i+k]
        if hamming_distance(window, pattern) <= max_mismatches:
            return i
    return -1

def extract_sequence(seq: str, up_seq: str, down_seq: str, max_mismatches: int) -> str:
    """
    Extract substring from seq between approximate matches of up_seq and down_seq.
    Parameters:
        -- seq: the sequence to search in
        -- up_seq: upstream sequence to match
        -- down_seq: downstream sequence to match
        -- max_mismatches: maximum number of mismatches allowed for both matches
    Returns:
        -- str: the extracted sequence or an error message if not found
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

def reverse_complement(seq: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.
    Parameters:
        -- seq: the DNA sequence to reverse complement
    Returns:
        -- str: the reverse complement of the sequence
    """
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def check_barcode(barcode_seq: str, barcode_temp: str, max_mismatches: int) -> bool:
    """
    Check if barcode_seq matches the barcode_temp allowing max_mismatches
    Parameters:
        -- barcode_seq: the sequence to check
        -- barcode_temp: the template sequence to match against
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- bool: True if matches within allowed mismatches, False otherwise
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
    Parameters:
        -- read1: tuple (header, sequence, quality) for read 1
        -- read2: tuple (header, sequence, quality) for read 2
        -- variant_up: upstream sequence for variant extraction
        -- variant_down: downstream sequence for variant extraction
        -- barcode_up: upstream sequence for barcode extraction
        -- barcode_down: downstream sequence for barcode extraction
        -- max_mismatches: maximum number of mismatches allowed for both extractions
    Returns:
        -- tuple: (variant_seq, barcode_seq) extracted sequences
    """
    seq1 = read1[1]
    seq2 = read2[1]

    variant_seq = extract_sequence(seq1, variant_up, variant_down, max_mismatches)
    barcode_seq = extract_sequence(seq2, barcode_up, barcode_down, max_mismatches)
    
    return (variant_seq, barcode_seq)

def batch_process_read_pairs(batch_reads, variant_up, variant_down, barcode_up, barcode_down, max_mismatches):
    """
    Process a batch of read pairs to extract variants and barcodes.
    Parameters:
        -- batch_reads: list of tuples, each containing (read1, read2)
        -- variant_up: upstream sequence for variant extraction
        -- variant_down: downstream sequence for variant extraction
        -- barcode_up: upstream sequence for barcode extraction
        -- barcode_down: downstream sequence for barcode extraction
        -- max_mismatches: maximum number of mismatches allowed for both extractions
    Returns:
        -- list of tuples: each tuple contains (variant_seq, barcode_seq) for processed read
    """
    results = []
    for r1, r2 in batch_reads:
        result = process_read_pair(r1, r2, variant_up, variant_down, barcode_up, barcode_down, max_mismatches)
        results.append(result)
    return results

def function_for_processpool(args):
    """
    Wrapper function for process pool as ProcessPoolExecutor expects a function rather than returned results.
    """
    return batch_process_read_pairs(*args)
    
def read_fastq_in_chunk(r1_path, r2_path, variant_up, variant_down, barcode_up, barcode_down, max_mismatches, chunk_size, processes):
    """
    Read paired-end FASTQ files in chunks and process each pair of reads in parallel
    Parameters:
        -- r1_path: path to read 1 FASTQ file
        -- r2_path: path to read 2 FASTQ file
        -- variant_up: upstream sequence for variant extraction
        -- variant_down: downstream sequence for variant extraction
        -- barcode_up: upstream sequence for barcode extraction
        -- barcode_down: downstream sequence for barcode extraction
        -- max_mismatches: maximum number of mismatches allowed for both extractions
        -- chunk_size: number of read pairs to process in each chunk
        -- processes: number of parallel processes to use
    Yields:
        -- list of tuples: each tuple contains (variant_seq, barcode_seq) for processed read
    """
    r1_handle = io.TextIOWrapper(pigz_open(r1_path).stdout) if r1_path.endswith(".gz") else open(r1_path)
    r2_handle = io.TextIOWrapper(pigz_open(r2_path).stdout) if r2_path.endswith(".gz") else open(r2_path)

    r1_iter = fastq_iter(r1_handle)
    r2_iter = fastq_iter(r2_handle)

    while True:
        r1_chunk = list()
        r2_chunk = list()
        try:
            for _ in range(chunk_size):
                r1_chunk.append(next(r1_iter))
                r2_chunk.append(next(r2_iter))
        except StopIteration:
            pass

        if not r1_chunk or not r2_chunk:
            break

        # Divide chunk into batches
        # Because process_read_pair is very fast, we can use a larger batch size
        # to make better use of CPU resources
        # if no batch, CPU is only 20% utilized
        batch_size = min(chunk_size, 5000)
        read_batches = [
            list(zip(r1_chunk[i:i+batch_size], r2_chunk[i:i+batch_size]))
            for i in range(0, len(r1_chunk), batch_size)
        ]

        args_list = [
            (batch, variant_up, variant_down, barcode_up, barcode_down, max_mismatches)
            for batch in read_batches
        ]

        with ProcessPoolExecutor(max_workers=processes) as executor:
            batch_results = list(executor.map(function_for_processpool, args_list))

        # Flatten list of lists
        results = [item for batch in batch_results for item in batch]
        yield results

        if len(r1_chunk) < chunk_size:
            break

    r1_handle.close()
    r2_handle.close()

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Extract variant and barcode from paired-end FASTQ files.", allow_abbrev = False)
    parser.add_argument("--read1",            type = str, required = True, help = "Read 1 FASTQ file")
    parser.add_argument("--read2",            type = str, required = True, help = "Read 2 FASTQ file")
    parser.add_argument("--variant_up",       type = str, required = True, help = "Upstream flank for variant sequence")
    parser.add_argument("--variant_down",     type = str, required = True, help = "Downstream flank for variant sequence")
    parser.add_argument("--variant_check",    action="store_true",         help = "Enable variant checking against length")
    parser.add_argument("--variant_len",      type = int,                  help = "Length of the variant sequence")
    parser.add_argument("--barcode_up",       type = str, required = True, help = "Sequence before barcode in read2")
    parser.add_argument("--barcode_down",     type = str, required = True, help = "Sequence after barcode in read2")
    parser.add_argument("--barcode_check",    action="store_true",         help = "Enable barcode checking against template")
    parser.add_argument("--barcode_temp",     type = str,                  help = "Template for barcode sequence (e.g., 'NNATNNNNATNNNNATNNNN')")
    parser.add_argument("--barcode_mismatch", type = int, default = 1,     help = "Number of mismatches allowed in barcode checking")
    parser.add_argument("--max_mismatches",   type = int, default = 2,     help = "Max mismatches allowed in up/down matches")
    parser.add_argument("--min_barcov",       type = int, default = 1,     help = "Minimum coverage for barcode-variant association")
    parser.add_argument("--chunk_size",       type = int, default = 5000,  help = "Chunk size for processing reads")
    parser.add_argument("--threads",          type = int, default = 12,    help = "Number of threads")
    parser.add_argument("--output_file",      type = str, default = "variant_barcode_association.txt", help = "Output file")
    parser.add_argument("--log_file",         type = str, default = "variant_barcode_association.log", help = "Log file")
    
    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.variant_check and args.variant_len is None:
        parser.error("--variant_len is required when --variant_check is set")

    if args.barcode_check and args.barcode_temp is None:
        parser.error("--barcode_temp is required when --barcode_check is set")

    print(f"Detecting variant and barcode associations, please wait...", flush=True)
    all_results = []
    for i, result_chunk in enumerate(read_fastq_in_chunk(args.read1, args.read2, 
                                                         args.variant_up.upper(), args.variant_down.upper(), 
                                                         args.barcode_up.upper(), args.barcode_down.upper(), 
                                                         args.max_mismatches, args.chunk_size, args.threads)):
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Processed chunk {i+1} with {len(result_chunk)} read pairs", flush=True)
        all_results.extend(result_chunk)
    print(f"Finished processing all the reads.", flush=True)

    count_varup_notfound = 0
    count_vardown_notfound = 0
    count_barup_notfound = 0
    count_bardown_notfound = 0
    if args.variant_check:
        count_variant_length = 0
    if args.barcode_check:
        count_barcode_pattern = 0
    count_effective_reads = 0
    all_results_filtered = []
    
    for (variant_seq, barcode_seq) in all_results:
        if variant_seq == "upstream not found":
            count_varup_notfound += 1
            continue

        if variant_seq == "downstream not found":
            count_vardown_notfound += 1
            continue

        if barcode_seq == "upstream not found":
            count_barup_notfound += 1
            continue

        if barcode_seq == "downstream not found":
            count_bardown_notfound += 1
            continue

        if args.variant_check and len(variant_seq) != args.variant_len:
            count_variant_length += 1
            continue

        if args.barcode_check and not check_barcode(reverse_complement(barcode_seq), args.barcode_temp.upper(), args.barcode_mismatch):
            count_barcode_pattern += 1
            continue

        count_effective_reads += 1
        all_results_filtered.append((variant_seq, reverse_complement(barcode_seq)))

    with open(args.log_file, "w") as f:
        f.write(f"Total reads processed: {len(all_results)}\n")
        f.write(f"Total reads with variant upstream not found: {count_varup_notfound}\n")
        f.write(f"Total reads with variant downstream not found: {count_vardown_notfound}\n")
        f.write(f"Total reads with barcode upstream not found: {count_barup_notfound}\n")
        f.write(f"Total reads with barcode downstream not found: {count_bardown_notfound}\n")
        if args.variant_check:
            f.write(f"Total reads with variant length inconsistent: {count_variant_length}\n")
        if args.barcode_check:
            f.write(f"Total reads with barcode pattern mismatch: {count_barcode_pattern}\n")
        f.write(f"Total effective reads: {count_effective_reads}\n")

    variant_barcode_association = Counter(all_results_filtered)
    with open(args.output_file, "w") as f:
        f.write("Barcode\tVariant\tCount\n")
        for (variant, barcode), count in variant_barcode_association.items():
            f.write(f"{barcode}\t{variant}\t{count}\n")
