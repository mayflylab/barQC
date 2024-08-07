#!/usr/bin/env python3

"""
Author: Maria Rossello
Date created: August 2024

Description:
    This script processes Split-seq reads to extract, correct, and tag barcodes. 
    The reads must be clean of Illumina adapters. Additionally, it generates 
    statistical reports and visualizations of the barcode distribution per Split-seq plate.

Usage:
    barQC.py -c <clean_reads.bam> -o <output.bam> [optional arguments]

Arguments:
    -c, --clean_reads : Path to the clean reads BAM file (required)
    -o, --output_bam  : Path to the output BAM file (required)
    -b, --bc_dir      : Directory with expected barcode files (default: current directory)
    -q, --qval        : Quality threshold (default: 10)
    -t, --threads     : Number of threads for parallel processing (default: 20)
    -v, --verbose     : Enable verbose logging

Output Files:
    <output.bam>              : BAM file with tagged barcodes.
    <output_base>.log         : Log file with processing information.
    <output_base>_debug.log   : Debug log file (if verbose mode is enabled).
    <output_base>_stats.log   : Statistics log file with processing metrics.
    <output_base>_heatmap_<barcode>.png : Heatmap visualizations of barcode distribution per Split-seq plate.
"""


#########################################################################################################
# PROGRAM ARGUMENTS
#########################################################################################################

#--------------------------------------------------------------------------------------------------------
# Load libraries to use in teh script
#--------------------------------------------------------------------------------------------------------

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
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import subprocess

# Suppress warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

#--------------------------------------------------------------------------------------------------------
# Argument parser setup
#--------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(prog='barQC.py',
                                 description='Process BAM files to extract, correct, and tag barcodes for Split-seq reads.',
                                 formatter_class=argparse.RawTextHelpFormatter)

# Group required arguments
required = parser.add_argument_group('required arguments')

required.add_argument('-c', '--clean_reads',
                      help='Path to the clean reads BAM file. Reads mus be without adapters.',
                      type=str,
                      required=True)

required.add_argument('-o', '--output_bam',
                      help='Path to the output BAM file.',
                      type=str,
                      required=True)

# Group optional arguments
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-b', '--bc_dir', 
                      help='Directory where the expected barcodes CSV files and linker fasta are stored. '
                           'Default to barcode_list.',
                      type=str, 
                      default="./barcode_list")

optional.add_argument('-q', '--qval', 
                      help='Quality threshold for barcode evaluation.',
                      type=int, 
                      default=10)

optional.add_argument('-t', '--threads',
                      help='Number of threads to use for parallel processing. Default is 20.',
                      type=int,
                      default=20)

optional.add_argument('-v', '--verbose', 
                      help='Enable verbose logging.',
                      action='store_true')

# Parse the arguments
args = parser.parse_args()

#--------------------------------------------------------------------------------------------------------
# Check all the requirements are provided
#--------------------------------------------------------------------------------------------------------

# Check if all required arguments are provided and validate paths
required_args_path = {'clean_reads': args.clean_reads, 'bc_dir': args.bc_dir}
missing_args_path = [arg for arg, path in required_args_path.items() if path is None or not os.path.isfile(path) and arg != 'bc_dir']
if missing_args_path:
    parser.print_help()
    print(f"\nError: Missing or invalid arguments or incorrect path: {', '.join(missing_args_path)}\n")
    exit(1)

# Check if bc_dir contains necessary files
expected_files = ['expected_barcodes_1.csv', 'expected_barcodes_2.csv', 'expected_barcodes_3.csv', 'invariable-linker-sequences.fasta']
missing_files = [file for file in expected_files if not os.path.isfile(os.path.join(args.bc_dir, file))]
if missing_files:
    parser.print_help()
    print(f"\nError: The following files are missing in the specified bc_dir ({args.bc_dir}): {', '.join(missing_files)}\n")
    exit(1)

# Check if all required arguments are provided
required_args = ['output_bam']
missing_args = [arg for arg in required_args if getattr(args, arg) is None]
if missing_args:
    parser.print_help()
    print(f"\nError: Missing required argument: {', '.join(missing_args)}\n")
    exit(1)


#--------------------------------------------------------------------------------------------------------
# Set up logging configurations
#--------------------------------------------------------------------------------------------------------

# Determine the output directory and base name from the output BAM file path
output_dir = os.path.dirname(args.output_bam)
output_base = os.path.basename(args.output_bam).split('.')[0]
log_path = os.path.join(output_dir, f"{output_base}.log")

# Configure the main log to log general information
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[
                        logging.FileHandler(log_path, mode='w') 
                    ])

# Configure the debug log if --verbose is provided
if args.verbose:
    debug_log_path = os.path.join(output_dir, f"{output_base}_debug.log")
    debug_logger = logging.getLogger('debug_logger')
    debug_logger.setLevel(logging.DEBUG)
    debug_handler = logging.FileHandler(debug_log_path, mode='w')
    debug_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    debug_logger.addHandler(debug_handler)
    debug_logger.propagate = False 
else:
    debug_logger = None

# Configure the logger for final statistics
stats_log_path = os.path.join(output_dir, f"{output_base}_stats.log")
stats_logger = logging.getLogger('stats_logger')
stats_logger.setLevel(logging.INFO)
stats_handler = logging.FileHandler(stats_log_path, mode='w')
stats_handler.setFormatter(logging.Formatter('%(message)s'))  # Only log the message
stats_logger.addHandler(stats_handler)
stats_logger.propagate = False


#########################################################################################################
# FUNCTIONS
#########################################################################################################

#--------------------------------------------------------------------------------------------------------
# Auxiliary functions to map the barcodes
#--------------------------------------------------------------------------------------------------------

def dataframe_to_fastq(bam_as_df):

    """
    Converts a DataFrame containing BAM data to FASTQ format.

    Args:
        bam_as_df (pd.DataFrame): DataFrame containing BAM data with columns 'query_name', 'seq', and 'qual'.

    Returns:
        str: FASTQ formatted string.
    """
   
    fastq_entries = []
    for _, row in bam_as_df.iterrows():
        query_name = row['query_name']
        seq = row['seq']
        qual = ''.join(chr(q + 33) for q in row['qual'])
        
        # Construct FASTQ entry
        fastq_entry = f"@{query_name}\n{seq}\n+\n{qual}\n"
        fastq_entries.append(fastq_entry)
    
    return ''.join(fastq_entries)

def trim_seq(fastq_data, debug_logger=None):

    """
    Trims sequences in FASTQ data to 100 bp using cutadapt.

    Args:
        fastq_data (str): FASTQ formatted string.
        debug_logger (logging.Logger, optional): Logger for debugging messages. Default is None.

    Returns:
        str: Trimmed FASTQ formatted string, or None if an error occurs.
    """

    try:
        process = subprocess.Popen(
            ['cutadapt', '-l', '100', '-'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        trimmed_fastq, cutadapt_stderr = process.communicate(input=fastq_data)
        
        if cutadapt_stderr and debug_logger:
            debug_logger.debug(f"Cutadapt:\n{cutadapt_stderr}")
        
        return trimmed_fastq
    
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred while running cutadapt: {e}")
        logging.error(f"Error message:\n{e.stderr}")
        return None
   
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
        tuple: Corrected barcode and a list of error codes, if any.
    """
    non_correct_bc = []
    
    if barcode in expected_barcode.values:
        # Perfect Barcode
        return barcode, non_correct_bc
    
    else:
        # Barcode need correction
        distances = [compute_hamming_distance(barcode, expected_barcode) for expected_barcode in expected_barcode.values]
        if None in distances:  
            # Different barcode length
            non_correct_bc.append('Err1')
            return None, non_correct_bc
        
        matching_indices = [i for i, distance in enumerate(distances) if distance == 1]
        if len(matching_indices) == 1:
            return expected_barcode.values[matching_indices[0]], non_correct_bc
        else:
            # It's impossible to correct the barcode
            non_correct_bc.append('Err2')
            return None, non_correct_bc


#--------------------------------------------------------------------------------------------------------
# Functions to open BAM file
#--------------------------------------------------------------------------------------------------------

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
                    if read.is_read2:
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
# Functions to obtain the align barcodes
#--------------------------------------------------------------------------------------------------------

def run_bbmap(bam_as_df, reference, cores):

    """
    Runs bbmap.sh to align sequences and find barcodes using in-memory FASTQ data.

    Args:
        bam_as_df (pd.DataFrame): DataFrame containing BAM data.
        reference (str): Path to the reference genome.
        cores (int): Number of threads to use for parallel processing.

    Returns:
        bytes: Mapped BAM data in SAM format.
    """
   
    logging.info("Obtaining FASTQ from BAM.")
    fastq = dataframe_to_fastq(bam_as_df)

    logging.info("Cutting sequences to 100 bp")
    trimmed_fastq = trim_seq(fastq)

    if trimmed_fastq is None:
        logging.error("Error: Timing did not run successfully.")
        sys.exit(1) 

    logging.info("Aligning sequence to find barcodes")

    try:
        process = subprocess.Popen(
            [
                'bbmap.sh',
                f'ref={reference}',
                'in=stdin.fastq',
                'out=stdout.sam',
                'int=f',
                'k=8',
                'minid=40',
                f'threads={cores}'
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=False  # Important for binary output
        )
        logging.debug("Communicating with bbmap.sh process.")
        bbmap_stdout, bbmap_stderr = process.communicate(input=trimmed_fastq.encode())

        logging.info("Sequence aligned to find barcodes")
        logging.debug(f"bbmap.sh stdout length: {len(bbmap_stdout)} bytes")
        logging.debug(f"bbmap.sh stderr length: {len(bbmap_stderr)} bytes")

        if bbmap_stderr:
            logging.debug(f"BBMap Info:\n{bbmap_stderr.decode()}")

        return bbmap_stdout
    
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred while running bbmap.sh: {e}")
        logging.error(f"Error message:\n{e.stderr.decode()}")
        sys.exit(1)

def parse_and_compute_coordinates(mapped_bam):

    """
    Parses the aligned BAM data to extract read names and CIGAR strings, 
    then computes the barcode coordinates based on CIGAR patterns.

    Args:
        mapped_bam (bytes): In-memory BAM data in SAM format.

    Returns:
        dict: Dictionary containing the computed barcode coordinates indexed by query names.
    """

    logging.info(f"Starting parsing barcodes to compute coordinates")

    records = []

    # Get CIGAR per sequence
    if isinstance(mapped_bam, bytes):
        mapped_bam = mapped_bam.decode('utf-8')
    lines = mapped_bam.split('\n')

    for line in lines:
        if not line.strip():
            # Skip empty lines
            continue

        if line.startswith('@'):
            # Skip header lines
            continue

        fields = line.split('\t')
        cigar = fields[5]

        if cigar and re.match(r'^[0-9]+S.*[^0-9][7-9]I.*', cigar):
            records.append({
                'query_name': fields[0].split(' ')[0],  
                'cigarstring': cigar
                })

    map_df = pd.DataFrame(records)

    # Get barcodes
    map_df = compute_initial_coordinates(map_df, 'BC1', r'([0-9]+S$)', 0)
    map_df = compute_initial_coordinates(map_df, 'BC2', r'([^0-9][7-9]I.*)', 0)
    map_df = compute_initial_coordinates(map_df, 'BC3', r'(?![0-9]).*$', -8)

    map_df = map_df[(map_df['BC1'] >= 17) & (map_df['BC2'] >= 9) & (map_df['BC3'] >= 1)]

    return map_df.set_index('query_name').to_dict('index')
                        

#--------------------------------------------------------------------------------------------------------
# Process barcodes
#--------------------------------------------------------------------------------------------------------

def process_barcodes(clean_df, map_dict, barcode_files_dir, qthreshold):

    """
    Processes the barcodes by extracting, evaluating quality, and correcting against expected barcodes.

    Args:
        clean_df (pd.DataFrame): DataFrame containing the clean reads data.
        map_dict (dict): Dictionary with mapping information.
        barcode_files_dir (str): Directory containing expected barcode files.
        qthreshold (int): Quality threshold for barcode evaluation.

    Returns:
        tuple: 
            pd.DataFrame: DataFrame containing the processed and corrected barcodes.
            pd.DataFrame: DataFrame with barcode statistics.
    """    
    
    total_rows = len(clean_df)
    barcode_ls = []
    
    try:
        bc1_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_1.csv"))
    except FileNotFoundError:
        logging.error("File 'expected_barcodes_1.csv' not found in the directory '%s'.", barcode_files_dir)
        return pd.DataFrame(), pd.DataFrame()

    try:
        bc2_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_2.csv"))
    except FileNotFoundError:
        logging.error("File 'expected_barcodes_2.csv' not found in the directory '%s'.", barcode_files_dir)
        return pd.DataFrame(), pd.DataFrame()

    try:
        bc3_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_3.csv"))
    except FileNotFoundError:
        logging.error("File 'expected_barcodes_3.csv' not found in the directory '%s'.", barcode_files_dir)
        return pd.DataFrame(), pd.DataFrame()

    # Initialize a DataFrame to store barcode counts
    barcode_stats = pd.DataFrame(columns=['WellPosition', 'Name', 'Barcode', 'Count'])

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
            def correct_and_append(barcode, qbarcode, expected_data, qname, bc_type):
                if evaluate_barcode_quality(qbarcode, qthreshold):
                    corrected, errors = check_and_correct_barcode(barcode, expected_data['Barcode'])
                    if debug_logger:
                        for error in errors:
                            if error == 'Err1':
                                debug_logger.debug(f"Wrong barcode length in {qname} for {bc_type}")
                            elif error == 'Err2':
                                debug_logger.debug(f"Impossible to correct the barcode in {qname} for {bc_type}")
                    return corrected, errors
                return None, []

            corrected_bc1, errors_bc1 = correct_and_append(bc1, qbc1, bc1_data, qname, 'Barcode1')
            corrected_bc2, errors_bc2 = correct_and_append(bc2, qbc2, bc2_data, qname, 'Barcode2')
            corrected_bc3, errors_bc3 = correct_and_append(bc3, qbc3, bc3_data, qname, 'Barcode3')

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

                # Update barcode statistics
                for corrected_bc, bc_data in zip([corrected_bc1, corrected_bc2, corrected_bc3], [bc1_data, bc2_data, bc3_data]):
                    if corrected_bc is not None:
                        well_position = bc_data.loc[bc_data['Barcode'] == corrected_bc, 'WellPosition'].values[0]
                        name = bc_data.loc[bc_data['Barcode'] == corrected_bc, 'Name'].values[0]
                        existing_row = barcode_stats[(barcode_stats['Barcode'] == corrected_bc) & (barcode_stats['WellPosition'] == well_position)]
                        if not existing_row.empty:
                            barcode_stats   .loc[existing_row.index, 'Count'] += 1
                        else:
                            new_row = pd.DataFrame([{'WellPosition': well_position, 'Name': name, 'Barcode': corrected_bc, 'Count': 1}])
                            barcode_stats = pd.concat([barcode_stats, new_row], ignore_index=True)

    return pd.DataFrame(barcode_ls), barcode_stats


#--------------------------------------------------------------------------------------------------------
# Make graphs and stats
#--------------------------------------------------------------------------------------------------------

def parse_well_position(df):

    """
    Parses the WellPosition column into separate Row and Column for plotting.

    Args:
        df (pd.DataFrame): DataFrame containing the barcode statistics.

    Returns:
        pd.DataFrame: Parsed DataFrame with separate Row and Column columns.
    """

    df['Row'] = df['WellPosition'].str[0]
    df['Column'] = df['WellPosition'].str[1:].astype(int)
    return df

def create_heatmap(data, title, filename):

    """
    Creates a heatmap for the given data and saves it to a file.

    Args:
        data (pd.DataFrame): DataFrame containing the parsed barcode statistics.
        title (str): Title for the heatmap.
        filename (str): Filename to save the heatmap.

    Returns:
        None
    """

    # Check if data is empty before pivoting
    if data.empty:
        logging.warning(f"No data available to create heatmap for {title}.")
        return

    # Pivot the data for heatmap and call infer_objects to prevent future warnings
    heatmap_data = data.pivot_table(index="Row", columns="Column", values="Count", fill_value=0).infer_objects(copy=False)
    
    # Sort the index in reverse alphabetical order
    heatmap_data = heatmap_data.loc[sorted(heatmap_data.index, reverse=True)]

    # Define custom colormap
    indigo_cmap = LinearSegmentedColormap.from_list('white_to_indigo', ['white', '#002395'], N=1000)
    
    # Create the heatmap
    plt.figure(figsize=(12, 10))  # Adjusted width for wider graph
    ax = sns.heatmap(heatmap_data, annot=True, fmt=".0f", cmap=indigo_cmap, cbar=False, 
                     xticklabels=True, yticklabels=True, linewidths=.5, linecolor='black')
    plt.title(title)
    plt.gca().invert_yaxis()
   
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.tick_params(left=False, top=False)
    
    # Save the heatmap to a file
    plt.savefig(filename)
    plt.close()


#--------------------------------------------------------------------------------------------------------
# Tag final bam
#--------------------------------------------------------------------------------------------------------

def tag_barcode(clean_reads, output_bam, corrected_barcode_df, num_threads, chunk_size=10000):
    
    """
    Tags barcodes in the clean reads BAM file.

    Args:
        clean_reads (str): Path to the clean reads BAM file.
        output_bam (str): Path to the output BAM file.
        corrected_barcode_df (pd.DataFrame): DataFrame containing the corrected barcodes.
        num_threads (int): Number of threads to use for parallel processing.
        chunk_size (int, optional): Number of reads to process in each chunk. Default is 10000.

    Returns:
        None
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


#########################################################################################################
# PROCESS DATA
#########################################################################################################

def main():
    # Some variables
    clean_reads = args.clean_reads
    num_threads = args.threads
    barcode_files_dir = args.bc_dir
    qthreshold = args.qval
    output_bam = args.output_bam
    reference = os.path.join(barcode_files_dir, "invariable-linker-sequences.fasta")

    logging.info("Starting barcode processing pipeline.")
    
    # First executor for initial parsing tasks
    with ThreadPoolExecutor(max_workers=num_threads) as initial_executor:
        future_clean_df = [initial_executor.submit(parse_clean_bam_file, clean_reads)]

        clean_df = None

        for future in as_completed(future_clean_df):
            try:
                clean_df = future.result()
                logging.info("Parsed Clean Reads BAM to DataFrame")
                gc.collect()
            except Exception as e:
                logging.error(f"Failed to parse clean BAM file: {e}")
                sys.exit(1)

    if clean_df is not None:
        # Use ThreadPoolExecutor for running bbmap and parsing coordinates
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            future_bbmap = executor.submit(run_bbmap, clean_df, reference, num_threads)

            try:
                mapped_bam_data = future_bbmap.result()
                logging.info("Barcode aligment to invariant sequence completed successfully.")

                # Submit the parse_and_compute_coordinates task
                future_map_dict = executor.submit(parse_and_compute_coordinates, mapped_bam_data)

                map_dict = future_map_dict.result()
                logging.info("Coordenated computed")
                gc.collect()

            except Exception as e:
                logging.error(f"Failed to align barcodes or parse coordinates: {e}")
                sys.exit(1)

    original_reads = len(clean_df) if clean_df is not None else 0
    tagged_reads = 0
    efficiency = 0.0
    num_cells = 0

    if map_dict and not clean_df.empty:
        # Second executor for barcode processing and tagging tasks
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            future_results = executor.submit(process_barcodes, clean_df, map_dict, barcode_files_dir=barcode_files_dir, qthreshold=qthreshold)

            try:
                barcode_df, stats_df = future_results.result()
                gc.collect()

                if not barcode_df.empty:
                    tagged_reads = len(barcode_df)
                    efficiency = (tagged_reads / original_reads) * 100 if original_reads > 0 else 0.0

                    unique_cells = barcode_df.drop_duplicates(subset=['BC1', 'BC2', 'BC3', 'UMI'])
                    num_cells = len(unique_cells)

                    parsed_df = parse_well_position(stats_df)
                    for barcode in parsed_df['Name'].unique():
                        barcode_data = parsed_df[parsed_df['Name'] == barcode]
                        heatmap_filename = os.path.join(output_dir, f"{output_base}_heatmap_{barcode}.png")
                        logging.info(f"Creating heatmap for {barcode}, filename: {heatmap_filename}")
                        create_heatmap(barcode_data, f"Heatmap for {barcode}", heatmap_filename)

                    tag_barcode(clean_reads, output_bam, barcode_df, num_threads=num_threads)
                    logging.info("Barcodes tagged to the final bam")
                    gc.collect()
                else:
                    logging.error("Processed barcode DataFrame is empty. Skipping tagging step.")
            except Exception as e:
                logging.error(f"Failed to process barcodes: {e}")
                sys.exit(1)

    # Print final statistics
    stats_logger.info(f"Number of original reads: {original_reads}")
    stats_logger.info(f"Number of tagged reads: {tagged_reads}")
    stats_logger.info(f"Efficiency: {efficiency:.2f}%")
    stats_logger.info(f"Number of unique cells: {num_cells}")

    # Print stats_df to stats_logger
    if not stats_df.empty:
        stats_df = stats_df.sort_values(by=['Name', 'WellPosition'])
        stats_df_dropped = stats_df.drop(columns=['Column', 'Row'])
        stats_logger.info("------\nBarcode Statistics:")
        stats_logger.info(stats_df_dropped.to_string(index=False))

    logging.info("Barcode processing pipeline completed.")
    logging.info("Have fun analyzing! (>ᴗ•)❀")

    # Ensure all logging is outputted
    for handler in logging.root.handlers[:]:
        handler.flush()
        handler.close()

if __name__ == "__main__":
    main()