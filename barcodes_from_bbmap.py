#!/usr/bin/env python3

"""
barcodes_from_bbmap_48plate.py [change?]

Author: Maria Rossello
Date created: May 2024

Description:
    This script processes mapped and cleaned BAM files to extract, correct, and tag barcodes. [+++]

Usage:
    python barcodes_from_bbmap_48plate.py -m <mapped_reads.bam> -c <clean_reads.bam> -o <output.bam> [-b <barcode_dir>] [-q <quality_threshold>] [-t <threads>]

Arguments:
    -m, --mapped_reads : Path to the mapped reads BAM file (output of bbmap).
    -c, --clean_reads  : Path to the clean reads BAM file (output of cutadapt).
    -o, --output_bam   : Path to the output BAM file.
    -b, --bc_dir       : Directory where the expected barcode files are stored.
    -q, --qval         : Quality threshold for evaluating barcode quality (default: 10).
    -t, --threads      : Number of threads to use for parallel processing (default: 20).

Requirements:
    - Python 3.9 or later
    - pandas
    - pysam
    - concurrent.futures
"""


#########################################################################################################
# PROGRAM ARGUMENTS
#########################################################################################################

import argparse
import os
import pandas as pd
import warnings
import pysam
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
from tqdm import tqdm
import numpy as np
import gc
import sys
import os
from queue import Queue

# Suppress warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Argument parser setup
parser = argparse.ArgumentParser(prog='barcodes_from_bbmap_48plate.py',
                                 description='Extract, correct, and tag barcodes from BAM files.',
                                 formatter_class=argparse.RawTextHelpFormatter)

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-m', '--mapped_reads',
                      help='Path to the mapped reads BAM file.',
                      type=str,
                      required=True)

required.add_argument('-c', '--clean_reads',
                      help='Path to the clean reads BAM file.',
                      type=str,
                      required=True)

required.add_argument('-o', '--output_bam',
                      help='Path to the output BAM file.',
                      type=str,
                      required=True)

optional.add_argument('-b', '--bc_dir', 
                      help='Directory where the expected barcode files are stored. Defaults to the directory this script is in.',
                      type=str, 
                      default=".")

optional.add_argument('-q', '--qval', 
                      help='Quality threshold',
                      type=int, 
                      default=10)

optional.add_argument('-t', '--threads',
                      help='Number of threads to use for parallel processing.',
                      type=int,
                      default=20)

args = parser.parse_args()


# Check if all required arguments are provided and validate paths
required_args_path = {'mapped_reads': args.mapped_reads, 'clean_reads': args.clean_reads}
missing_args_path = [arg for arg, path in required_args_path.items() if path is None or not os.path.isfile(path)]
if missing_args_path:
    parser.print_help()
    print(f"\nError: Missing or invalid arguments or incorrect path: {', '.join(missing_args_path)}\n")
    exit(1)

# Check if all required arguments are provided
required_args = ['output_bam']
missing_args = [arg for arg in required_args if getattr(args, arg) is None]
if missing_args:
    parser.print_help()
    print(f"\nError: Missing required argument: {', '.join(missing_args)}\n")
    exit(1)


logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[
                        logging.StreamHandler(sys.stdout)
                    ])


#########################################################################################################
# FUNCTIONS
#########################################################################################################

#--------------------------------------------------------------------------------------------------------
# Functions to map the barcodes from the mapping
#--------------------------------------------------------------------------------------------------------


def compute_initial_coordinates(map_df, barcode_col, pattern, offset):

    """
    Computes initial coordinates for barcodes based on a given pattern.

    Args:
        map_df (pd.DataFrame): DataFrame with mapped data.
        barcode_col (str): Column name to store the barcode coordinates.
        pattern (str): Regex pattern to identify barcode regions.
        offset (int): Offset to adjust the coordinates.

    Returns:
        pd.DataFrame: Updated DataFrame with computed barcode coordinates.
    """

    map_df[barcode_col] = map_df['cigarstring'].str.replace(pattern, '', regex=True)
    map_df[barcode_col] = map_df[barcode_col].str.replace(r'[0-9]+D', '', regex=True)
    map_df[barcode_col] = map_df[barcode_col].str.findall(r'\d+').apply(np.array, dtype=int).apply(np.sum) + offset
    map_df[barcode_col] = map_df[barcode_col].abs()
    
    return map_df


#--------------------------------------------------------------------------------------------------------
# Auxiliary functions to handle barcodes
#--------------------------------------------------------------------------------------------------------


def evaluate_barcode_quality(quality_scores, threshold):

    """
    Evaluates the quality of a barcode based on its quality scores.

    Args:
        quality_scores (list): List of quality scores.
        threshold (int): Quality threshold.

    Returns:
        bool: True if the average quality meets or exceeds the threshold, False otherwise.
    """

    return len(quality_scores) > 0 and np.mean(np.array(quality_scores)) >= threshold

def compute_hamming_distance(mybc, expectedbc):
    """
    Calculates the Hamming distance between two strings.

    Args:
        mybc (str): Barcode to be evaluated.
        expectedbc (str): Expected barcode.

    Returns:
        int: Hamming distance, or None if the barcodes have different lengths.
    """
    if len(mybc) != len(expectedbc):
        return None

    mybc_array = np.frombuffer(mybc.encode(), dtype=np.uint8)
    expectedbc_array = np.frombuffer(expectedbc.encode(), dtype=np.uint8)

    return np.sum(mybc_array != expectedbc_array)

def check_and_correct_barcode(barcode, expected_barcode):
    """
    Checks and corrects a barcode against a list of expected barcodes.

    Args:
        barcode (str): Barcode to be checked.
        expected_barcode (pd.Series): Series of expected barcodes.

    Returns:
        tuple: Corrected barcode and error code, if any.
    """
    if barcode in expected_barcode.values:
        #Perfect Barcode
        return barcode, None
    
    else:
        #Barcode need correction
        distances = [compute_hamming_distance(barcode, expected_barcode) for expected_barcode in expected_barcode.values]
        if None in distances:  
            # Different barcode length
            return None, 'Err1'
        
        matching_indices = [i for i, distance in enumerate(distances) if distance == 1]
        if len(matching_indices) == 1:
            return expected_barcode.values[matching_indices[0]], None
        else:
            # It's impossible to correct the barcode
            return None, 'Err2'


#--------------------------------------------------------------------------------------------------------
# Functions to open BAM files
#--------------------------------------------------------------------------------------------------------


def parse_and_compute_coordinates(mapped_bam_file_path):

    """
    Parses the alignment BAM file to extract read name and CIGAR, filtering out records that do not meet the minimal requirements. 
    Then it transforms the cigar information into coordinates to find the different barcode.

    Args:
        mapped_bam_file_path (str): Path to the mapped reads BAM file.

    Returns:
        pd.DataFrame: DataFrame containing the parsed and computed information.
    """

    logging.info(f"Opening Mapped Reads bam")
    records = []
    try:
        with pysam.AlignmentFile(mapped_bam_file_path, "rb") as bamfile:
            total_reads = bamfile.mapped + bamfile.unmapped
            with tqdm(total=total_reads, desc="Parsing mapped BAM file") as pbar:
                for read in bamfile.fetch(until_eof=True):
                    if read.flag == 137:
                        # Process only read2
                        
                        cigar = read.cigarstring
                        if cigar and re.match(r'^[0-9]+S.*[^0-9][7-9]I.*', cigar):
                            # The minimal accepted:
                            #   Some nucleotides before the invariant seq 
                            #   7-9 nt between invariant sequences
                            
                            records.append({
                                'query_name': read.query_name.split(' ')[0],
                                'cigarstring': cigar
                            })
                        pbar.update(1)

        map_df = pd.DataFrame(records)

        # Get barcodes
        map_df = compute_initial_coordinates(map_df, 'BC1', r'([0-9]+S$)', 0)
        map_df = compute_initial_coordinates(map_df, 'BC2', r'([^0-9][7-9]I.*)', 0)
        map_df = compute_initial_coordinates(map_df, 'BC3', r'(?![0-9]).*$', -8)

        map_df = map_df[(map_df['BC1'] >= 17) & (map_df['BC2'] >= 9) & (map_df['BC3'] >= 1)]

        return map_df.set_index('query_name').to_dict('index')
                        
    except Exception as e:
        logging.exception(f"Failed to parse BAM file {mapped_bam_file_path}: {e}")
        
        return {}

def parse_clean_bam_file(clean_bam_file_path):

    """
    Parses the clean BAM file to extract relevant information.

    Args:
        clean_bam_file_path (str): Path to the clean reads BAM file.

    Returns:
        pd.DataFrame: DataFrame containing the parsed information.
    """

    logging.info(f"Opening Cleaned Reads bam")

    records = []
    try:
        with pysam.AlignmentFile(clean_bam_file_path, "rb", check_sq=False) as bamfile:
            total_reads = bamfile.mapped + bamfile.unmapped
            with tqdm(total=total_reads, desc="Parsing clean BAM file") as pbar:
                for read in bamfile.fetch(until_eof=True):
                    if read.flag == 141:
                        # Process only read2
                            
                        records.append({
                            'query_name': read.query_name.split(' ')[0],
                            'seq': read.query_sequence,
                            'qual':read.query_qualities
                        })
                    pbar.update(1)    

        return pd.DataFrame(records)
    except Exception as e:
        logging.error(f"Failed to parse BAM file {clean_bam_file_path}: {e}")
        return pd.DataFrame()


#--------------------------------------------------------------------------------------------------------
# Process barcodes
#--------------------------------------------------------------------------------------------------------

# def process_barcodes(clean_df, map_dict, barcode_files_dir, qthreshold, chunk_size=100000):
#     """
#     Processes the barcodes by extracting, evaluating quality, and correcting against expected barcodes.
#     """
#     total_rows = len(clean_df)
#     barcode_ls = []

#     bc1_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_1.csv"))
#     bc2_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_2.csv"))
#     bc3_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_3.csv"))

#     def process_chunk(chunk):
#         local_barcode_ls = []
#         for entry in chunk.itertuples():
#             try:
#                 if entry.query_name in map_dict:
#                     qname = entry.query_name
#                     rmap = map_dict[qname]
#                     seq = np.array(list(entry.seq))
#                     q = np.array(entry.qual)

#                     bc1 = ''.join(seq[rmap['BC1']:rmap['BC1'] + 8])
#                     bc2 = ''.join(seq[rmap['BC2']:rmap['BC2'] + 8])
#                     bc3 = ''.join(seq[rmap['BC3']:rmap['BC3'] + 8])

#                     qbc1 = q[rmap['BC1']:rmap['BC1'] + 8]
#                     qbc2 = q[rmap['BC2']:rmap['BC2'] + 8]
#                     qbc3 = q[rmap['BC3']:rmap['BC3'] + 8]

#                     def correct_and_append(barcode, qbarcode, expected_data):
#                         if evaluate_barcode_quality(qbarcode, qthreshold):
#                             corrected, error = check_and_correct_barcode(barcode, expected_data['Barcode'])
#                             if error == 'Err1':
#                                 logging.error(f"Error: Barcode is not of the correct length in {qname}")
#                             return corrected
#                         return None

#                     corrected_bc1 = correct_and_append(bc1, qbc1, bc1_data)
#                     corrected_bc2 = correct_and_append(bc2, qbc2, bc2_data)
#                     corrected_bc3 = correct_and_append(bc3, qbc3, bc3_data)

#                     umi_start = rmap['BC3'] - 11 if rmap['BC3'] > 10 else 0
#                     umi_end = rmap['BC3'] if rmap['BC3'] > 1 else 1
#                     qumi = q[umi_start:umi_end]
#                     umi = ''.join(seq[umi_start:umi_end])

#                     if evaluate_barcode_quality(qumi, qthreshold):
#                         umi = umi
#                     else:
#                         umi = None

#                     if None not in [corrected_bc1, corrected_bc2, corrected_bc3]:
#                         local_barcode_ls.append({
#                             'query_name': qname,
#                             'BC1': corrected_bc1,
#                             'BC2': corrected_bc2,
#                             'BC3': corrected_bc3,
#                             'UMI': umi
#                         })
#             except Exception as e:
#                 logging.error(f"Error processing entry {entry}: {e}")
#         return local_barcode_ls

#     with ThreadPoolExecutor() as executor:
#         futures = []
#         for start in range(0, total_rows, chunk_size):
#             chunk = clean_df.iloc[start:start+chunk_size]
#             futures.append(executor.submit(process_chunk, chunk))

#         for future in as_completed(futures):
#             barcode_ls.extend(future.result())

#     return pd.DataFrame(barcode_ls)


def process_barcodes(clean_df, map_dict, barcode_files_dir, qthreshold):
    """
    Processes the barcodes by extracting, evaluating quality, and correcting against expected barcodes.

    Args:
        clean_df (pd.DataFrame): DataFrame containing the clean reads data.
        map_dict (dict): Dictionary with mapping information.
        barcode_files_dir (str): Directory containing expected barcode files.
        qthreshold (int): Quality threshold for barcode evaluation.

    Returns:
        pd.DataFrame: DataFrame containing the processed and corrected barcodes.
    """
    total_rows = len(clean_df)
    barcode_ls = []

    try:
        bc1_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_1.csv"))
    except FileNotFoundError:
        logging.error("File 'expected_barcodes_1.csv' not found in the directory '%s'.", barcode_files_dir)
        return bc1_data

    try:
        bc2_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_2.csv"))
    except FileNotFoundError:
        logging.error("File 'expected_barcodes_2.csv' not found in the directory '%s'.", barcode_files_dir)
        return bc2_data

    try:
        bc3_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_3.csv"))
    except FileNotFoundError:
        logging.error("File 'expected_barcodes_3.csv' not found in the directory '%s'.", barcode_files_dir)
        return bc3_data

    for entry in tqdm(clean_df.itertuples(), total=total_rows, desc="Processing and correcting barcodes"):

        if entry.query_name in map_dict:
        # Only process files that have a correct CIGAR

            qname = entry.query_name
            rmap = map_dict[qname]
            seq = np.array(list(entry.seq))
            q = np.array(entry.qual)

            # Get barcode sequence
            bc1 = ''.join(seq[rmap['BC1']:rmap['BC1'] + 10])
            bc2 = ''.join(seq[rmap['BC2']:rmap['BC2'] + 8])
            bc3 = ''.join(seq[rmap['BC3']:rmap['BC3'] + 8])

            # Get barcode quality score
            qbc1 = q[rmap['BC1']:rmap['BC1'] + 10]
            qbc2 = q[rmap['BC2']:rmap['BC2'] + 8]
            qbc3 = q[rmap['BC3']:rmap['BC3'] + 8]

            # Filter by quality and correct barcode
            def correct_and_append(barcode, qbarcode, expected_data):
                if evaluate_barcode_quality(qbarcode, qthreshold):
                    corrected, error = check_and_correct_barcode(barcode, expected_data['Barcode'])
                    if error == 'Err1':
                        logging.error(f"Error: Barcode is not of the correct length in {qname}")
                    return corrected
                return None

            corrected_bc1 = correct_and_append(bc1, qbc1, bc1_data)
            corrected_bc2 = correct_and_append(bc2, qbc2, bc2_data)
            corrected_bc3 = correct_and_append(bc3, qbc3, bc3_data)

            # Get UMI quality
            umi_start = rmap['BC3'] - 11 if rmap['BC3'] > 10 else 0
            umi_end = rmap['BC3'] if rmap['BC3'] > 1 else 1
            qumi = q[umi_start:umi_end]
            umi = ''.join(seq[umi_start:umi_end])

            # Get UMI if it has good quality
            if evaluate_barcode_quality(qumi, qthreshold):
                umi = umi
            else:
                umi = None

            if None not in [corrected_bc1, corrected_bc2, corrected_bc3]:
                barcode_ls.append({
                    'query_name': qname,
                    'BC1': corrected_bc1,
                    'BC2': corrected_bc2,
                    'BC3': corrected_bc3,
                    'UMI': umi
                })

    return pd.DataFrame(barcode_ls)


#--------------------------------------------------------------------------------------------------------
# Tag final bam
#--------------------------------------------------------------------------------------------------------

def tag_barcode(clean_reads, output_bam, corrected_barcode_df, num_threads, chunk_size=10000):
    """
    Tags barcodes in the clean reads BAM file.
    """

    def tag_barcodes_chunk(chunk, barcode_dict):
        output_reads = []
        for read in chunk:
            if read.flag == 77:
                qname = read.query_name.split(' ')[0]
                if qname in barcode_dict:
                    barcode_info = barcode_dict[qname]
                    xd, xe, xf, xm = barcode_info['BC1'], barcode_info['BC2'], barcode_info['BC3'], barcode_info['UMI']
                    xc = xd + xe + xf

                    # Add tags to the read
                    read.set_tag('XD', xd, value_type='Z')
                    read.set_tag('XE', xe, value_type='Z')
                    read.set_tag('XF', xf, value_type='Z')
                    read.set_tag('XC', xc, value_type='Z')
                    read.set_tag('XM', xm, value_type='Z')
                    output_reads.append(read)
        return output_reads

    def worker(bam_queue, result_queue, barcode_dict):
        while not bam_queue.empty():
            chunk = bam_queue.get()
            if chunk is None:
                break
            result_queue.put(tag_barcodes_chunk(chunk, barcode_dict))
            bam_queue.task_done()

    barcode_dict = corrected_barcode_df.set_index('query_name').to_dict('index')

    logging.info("Tagging barcodes in clean reads BAM file.")
    try:
        bam_queue = Queue()
        result_queue = Queue()

        with pysam.AlignmentFile(clean_reads, "rb", check_sq=False) as bamfile, \
                pysam.AlignmentFile(output_bam, "wb", header=bamfile.header) as outfile:
            
            total_reads = bamfile.mapped + bamfile.unmapped
            chunk = []
            for read in tqdm(bamfile.fetch(until_eof=True), total=total_reads, desc="Tagging barcodes"):
                chunk.append(read)
                if len(chunk) >= chunk_size:
                    bam_queue.put(chunk)
                    chunk = []
            if chunk:
                bam_queue.put(chunk)

            threads = []
            for _ in range(num_threads):
                thread = ThreadPoolExecutor(max_workers=1).submit(worker, bam_queue, result_queue, barcode_dict)
                threads.append(thread)
            
            bam_queue.join()
            
            while not result_queue.empty():
                tagged_reads = result_queue.get()
                for tagged_read in tagged_reads:
                    outfile.write(tagged_read)
                result_queue.task_done()
                
            for thread in threads:
                thread.result()  # Ensures all threads complete

    except Exception as e:
        logging.error(f"Failed to tag barcodes: {e}")


# def tag_barcode(clean_reads, output_bam, corrected_barcode_df, num_threads, chunk_size=10000):
#     """
#     Tags barcodes in the clean reads BAM file.

#     Args:
#         clean_reads (str): Path to the clean reads BAM file.
#         output_bam (str): Path to the output BAM file.
#         corrected_barcode_df (pd.DataFrame): DataFrame with corrected barcodes.
#         num_threads (int): Number of threads to use for parallel processing.
#         chunk_size (int): Size of the chunks for processing.

#     Returns:
#         None
#     """
#     def tag_barcodes_chunk(chunk, corrected_barcode_df):
#         output_reads = []
#         for read in chunk:
#             if read.flag == 77:
#                 qname = read.query_name.split(' ')[0]
#                 row = corrected_barcode_df[corrected_barcode_df['query_name'] == qname]

#                 if not row.empty:
#                     xd = row['BC1'].values[0]
#                     xe = row['BC2'].values[0]
#                     xf = row['BC3'].values[0]
#                     xc = xd + xe + xf
#                     xm = row['UMI'].values[0]

#                     # Add tags to the read
#                     read.set_tag('XD', xd, value_type='Z')
#                     read.set_tag('XE', xe, value_type='Z')
#                     read.set_tag('XF', xf, value_type='Z')
#                     read.set_tag('XC', xc, value_type='Z')
#                     read.set_tag('XM', xm, value_type='Z')
#                     output_reads.append(read)
#         return output_reads

#     logging.info("Tagging barcodes in clean reads BAM file.")
#     try:
#         with pysam.AlignmentFile(clean_reads, "rb", check_sq=False) as bamfile, \
#             pysam.AlignmentFile(output_bam, "wb", header=bamfile.header) as outfile:
        
#             total_reads = bamfile.mapped + bamfile.unmapped
#             chunk = []
#             futures = []
#             with ThreadPoolExecutor(max_workers=num_threads) as executor:
#                 for read in tqdm(bamfile.fetch(until_eof=True), total=total_reads, desc="Tagging barcodes"):
#                     chunk.append(read)
#                     if len(chunk) >= chunk_size:
#                         futures.append(executor.submit(tag_barcodes_chunk, chunk, corrected_barcode_df))
#                         chunk = []
#                 if chunk:
#                     futures.append(executor.submit(tag_barcodes_chunk, chunk, corrected_barcode_df))

#                 for future in as_completed(futures):
#                     for tagged_read in future.result():
#                         outfile.write(tagged_read)
#     except Exception as e:
#         logging.error(f"Failed to tag barcodes: {e}")


#########################################################################################################
# PROCESS DATA
#########################################################################################################

def main():

    # some variabes
    mapped_reads=args.mapped_reads
    clean_reads=args.clean_reads
    num_threads = args.threads
    barcode_files_dir=args.bc_dir
    qthreshold=args.qval
    output_bam=args.output_bam

    logging.info("Starting barcode processing pipeline.")
    
    # First executor for initial parsing tasks
    with ThreadPoolExecutor(max_workers=num_threads) as initial_executor:
        future_map_dict = initial_executor.submit(parse_and_compute_coordinates, mapped_reads)
        future_clean_df = initial_executor.submit(parse_clean_bam_file, clean_reads)

        map_dict = None
        clean_df = None

        for future in as_completed([future_map_dict, future_clean_df]):
            if future == future_map_dict:
                map_dict = future.result()
                logging.info("Parsed and computed coordinates from Mapped Reads BAM")
                gc.collect()
            elif future == future_clean_df:
                clean_df = future.result()
                logging.info("Parsed Clean Reads BAM to DataFrame")
                gc.collect()

    if map_dict and not clean_df.empty:
        # Second executor for barcode processing and tagging tasks
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            future_barcode_df = executor.submit(process_barcodes, clean_df, map_dict, barcode_files_dir=barcode_files_dir, qthreshold=qthreshold)

            barcode_df = future_barcode_df.result()
            gc.collect()
            
            if not barcode_df.empty:
                # Ensuring tag_barcode completes before shutting down the executor
                tag_barcode(clean_reads, output_bam, barcode_df, num_threads=num_threads)
                logging.info("Barcode processing pipeline completed successfully.")
                logging.info("Have fun analyzing! (>ᴗ•)❀")
                gc.collect()
            else:
                logging.error("Processed barcode DataFrame is empty. Skipping tagging step.")

    # Ensure all logging is outputted
    for handler in logging.root.handlers[:]:
        handler.flush()
        handler.close()

if __name__ == "__main__":
    main()
