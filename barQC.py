#!/usr/bin/env python3

"""
Author: Maria Rossello
Date created: August 2024

Description:
    This script processes Split-seq reads to extract, correct, and tag barcodes. 
    The reads must be clean of Illumina adapters. Additionally, it generates 
    statistical reports and visualizations of the barcode distribution per Split-seq plate.

Usage:
    barQC.py -f1 <read1.fastq> -f2 <read2.fastq> -o <output.bam> [optional arguments]

Arguments:
    -f1, --read1_fastq : Path to the fastq file of read1 (required)
    -f2, --read2_fastq : Path to the fastq file of read2 containing the barcode sequence (required)
    -o, --output_name   : Path to the output files (required)
    -b, --bc_dir       : Directory with expected barcode files (default: current directory)
    -q, --qval         : Quality threshold (default: 10)
    -t, --threads      : Number of threads for parallel processing (default: 20)
    -s, --stats        : Save stats in an output file
    -v, --verbose      : Enable verbose logging
    --skip_tagging     : Skip the tagging step and only generate statistics.

Output Files:
    <output_name>.bam                   : BAM file with tagged barcodes.
    <output_name>.log                   : Log file with processing information.
    <output_name>_heatmap_<barcode>.png : Heatmap visualizations of barcode distribution per Split-seq plate.
    <output_name>_stats.log             : Statistics log file with processing metrics (if stats mode is enabled).
    <output_name>_debug.log             : Debug log file (if verbose mode is enabled).
"""


#%% PROGRAM ARGUMENTS
#########################################################################################################

#%%% Load libraries to use in the script
#--------------------------------------------------------------------------------------------------------

import argparse
import os
import pandas as pd
import warnings
import pysam
import re
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
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
import psutil
import resource

# Suppress warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)


#%%% Argument parser setup
#--------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(prog='barQC.py',
                                 description='Process BAM files to extract, correct, and tag barcodes for Split-seq reads.',
                                 formatter_class=argparse.RawTextHelpFormatter)

# Group required arguments
required = parser.add_argument_group('required arguments')

required.add_argument('-f1', '--read1_fastq',
                      help='Path to the fastq file for read 1. This reads contains mRNA information' 
                           'Reads must be without adapters.',
                      type=str,
                      required=True)

required.add_argument('-f2', '--read2_fastq',
                      help='Path to the fastq file for read 2. This reads contains barcode information' 
                           'Reads must be without adapters.',
                      type=str,
                      required=True)

required.add_argument('-o', '--output_name',
                      help='Path to the output files. Add path and prefix.',
                      type=str,
                      required=True)

# Group optional arguments
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-b', '--bc_dir', 
                      help='Directory where the expected barcodes CSV files and linker fasta are stored.'
                           'Defaults to the directory this script is in.',
                      type=str, 
                      default=".")

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

optional.add_argument('-s', '--stats', 
                      help='Save stats in an output file.',
                      action='store_true')

optional.add_argument('--skip_tagging', 
                      help='Skip the tagging step and only generate statistics.', 
                      action='store_true')

optional.add_argument('--memory_limit', 
                      help='Limit memory usage (in GB). If not provided, the system will estimate a safe limit.', 
                      type=int, 
                      default=0)

# Parse the arguments
args = parser.parse_args()

#%%% Set up logging configurations
#--------------------------------------------------------------------------------------------------------

# Determine the output directory and base name from the output BAM file path
output_dir = os.path.dirname(args.output_name)
output_base = os.path.basename(args.output_name)
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
if args.stats:
    stats_log_path = os.path.join(output_dir, f"{output_base}_stats.log")
    stats_logger = logging.getLogger('stats_logger')
    stats_logger.setLevel(logging.INFO)
    stats_handler = logging.FileHandler(stats_log_path, mode='w')
    stats_handler.setFormatter(logging.Formatter('%(message)s'))  # Only log the message
    stats_logger.addHandler(stats_handler)
    stats_logger.propagate = False
else:
    stats_logger = None

#%%% Check all the requirements are provided
#--------------------------------------------------------------------------------------------------------

# Check if all required arguments are provided and validate paths
required_args_path = {'read1_fastq': args.read1_fastq, 'read2_fastq': args.read2_fastq,'bc_dir': args.bc_dir}
missing_args_path = [arg for arg, path in required_args_path.items() if (path is None or not os.path.isfile(path)) and arg != 'bc_dir']
if missing_args_path:
    parser.print_help()
    print(f"\n\nERROR: Missing or invalid arguments or incorrect path: {', '.join(missing_args_path)}\n")
    logging.error(f"Missing or invalid arguments or incorrect path: {', '.join(missing_args_path)}\n")
    exit(1)


# Check if bc_dir contains necessary files
expected_files = ['expected_barcodes_1.csv', 'expected_barcodes_2.csv', 'expected_barcodes_3.csv', 'invariable-linker-sequences.fasta']
missing_files = [file for file in expected_files if not os.path.isfile(os.path.join(args.bc_dir, file))]
if missing_files:
    parser.print_help()
    print(f"\n\nERROR: The following files are missing in the specified bc_dir ({args.bc_dir}): {', '.join(missing_files)}\n")
    logging.error(f"The following files are missing in the specified bc_dir ({args.bc_dir}): {', '.join(missing_files)}\n")
    exit(1)

# Check if all required arguments are provided
required_args = ['output_name']
missing_args = [arg for arg in required_args if getattr(args, arg) is None]
if missing_args:
    parser.print_help()
    print(f"\n\nERROR: Missing required argument: {', '.join(missing_args)}\n")
    logging.error(f"\n\nERROR: Missing required argument: {', '.join(missing_args)}\n")
    exit(1)

# Check if bbmap.sh is installed and accessible
try:
    result = subprocess.run(['bbmap.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Check the return code to confirm if bbmap.sh exists
    if result.returncode == 0:
        debug_logger.debug("bbmap.sh is installed and accessible.")
    else:
        logging.error(
            "bbmap.sh is not installed or not accessible."
            f"\nError: {result.stderr}"
            )
        
        print(f"\nError: bbmap.sh is not installed or not accessible. Exiting.\n")
        exit(1)

except FileNotFoundError:
    logging.error("bbmap.sh is not installed or not found in the system's PATH.")
    print(f"\nError: bbmap.sh is not installed or not found in the system's PATH. Exiting.\n")
    exit(1)

except Exception as e:
    logging.error(f"An unexpected error occurred while checking bbmap.sh: {e}")
    print(f"\nERROR: An unexpected error occurred while checking bbmap.sh. Exiting.\n")
    exit(1)


#%% FUNCTIONS
#########################################################################################################

#%%% Auxiliary functions calculate memory limitations
#--------------------------------------------------------------------------------------------------------

def estimate_memory_per_record_fq(file_path):

    file_size_bytes = os.path.getsize(file_path)
    with open(file_path, 'r') as f:
        num_lines = sum(1 for _ in f)
        num_records = num_lines // 4
    
    if num_records == 0:
        return 0
    
    memory_per_record_bytes = file_size_bytes / num_records
    memory_per_record_mb = memory_per_record_bytes / (1024 * 1024)
    
    return memory_per_record_mb

def calculate_available_memory(usr_memory_limit):

    if usr_memory_limit > 0:
        system_memory = psutil.virtual_memory().available // (1024 * 1024)
        available_memory = min(usr_memory_limit, system_memory)
    else:
        available_memory = psutil.virtual_memory().available // (1024 * 1024)

    return available_memory

def calculate_jvm_memory(available_memory, safety_margin=0.7):
    """
    Calculate JVM memory allocation for bbmap.sh based on available memory.
    
    Args:
        available_memory (int): Available memory in MB.
        safety_margin (float): Fraction of total memory to allocate to JVM (default is 70%).

    Returns:
        str: Memory allocation for JVM in the format required by BBMap (e.g., '4g' for 4 GB).
    """

    jvm_memory_mb = int(available_memory * safety_margin)
    
    jvm_memory_gb = jvm_memory_mb // 1024
    return f'{jvm_memory_gb}g'

def estimate_memory_per_record_df(df):

    num_records = len(df)
    total_size_bytes = df.memory_usage(deep=True).sum()
    memory_per_record = total_size_bytes / num_records if num_records > 0 else 0
    return memory_per_record / (1024 * 1024)

def calculate_chunk_size(available_memory, memory_per_record, safety_margin=0.7, threads=1):
    """
    Calculate the chunk size based on available memory and memory per record,
    reassessing before each process to adapt to real-time memory conditions.
    
    Args:
        available_memory (int): Available memory in MB.
        memory_per_record (float): Estimated memory usage per record in MB.
        safety_margin (float): Fraction of memory to actually use (default is 70%).
    
    Returns:
        int: Dynamically calculated chunk size.
    """
    usable_memory = (available_memory * safety_margin) / threads
    chunk_size = int(usable_memory / memory_per_record)
    return max(chunk_size, 1000)

def calculate_dynamic_process_count(available_memory, process_memory_overhead=500, safety_margin=0.7):
    """
    Calculate the number of processes based on available memory and an overhead estimate for each process.
    
    Args:
        available_memory (int): Available memory in MB.
        process_memory_overhead (int): Estimated memory overhead for each process in MB (default is 500MB).
        safety_margin (float): Fraction of memory to use (default is 70%).

    Returns:
        int: Maximum number of processes that can run without exhausting memory.
    """

    # Get current file descriptor limit
    
    # Calculate based on memory
    usable_memory = available_memory * safety_margin
    max_processes_by_memory = int(usable_memory / process_memory_overhead)

    # Calculate based on file descriptors
    soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
    max_processes_by_fd = soft_limit // 10

    return max(1, min(max_processes_by_memory, max_processes_by_fd))


#%%% Functions handel fasta files
#--------------------------------------------------------------------------------------------------------

def fastq_to_dataframe(fastq_file):
    query_names = []
    sequences = []
    qualities = []

    string_fastq = []

    try:
        # Calculate the total number of reads
        with pysam.FastxFile(fastq_file) as fastq:
            total_reads = sum(1 for _ in fastq)  # Count entries

        debug_logger.debug(
            f"Number of entries detected in the fastq file: {total_reads}"
        )

        # Parsing fastq with progress bar
        with pysam.FastxFile(fastq_file) as fastq:
            with tqdm(total=total_reads, desc="Parsing fastq for read 2") as pbar:
                for entry in fastq:
                    query_names.append(entry.name)
                    sequences.append(entry.sequence)
                    qualities.append(entry.quality)
                    pbar.update(1)

                    # Fastq to string for the 100 nucleotides
                    string_fastq.append(f"@{entry.name}\n{entry.sequence[:100]}\n+\n{entry.quality[:100].strip()}")

        # Create DataFrame
        df_fastq = pd.DataFrame({
            'query_name': query_names,
            'seq': sequences,
            'qual': qualities
        })

        debug_logger.debug(
            f"Number of entries converted to dataframe: {len(df_fastq)} "
            "(Fun: fastq_to_dataframe)"
        )

        debug_logger.debug(
            f"Number of entries converted to string (for bbmap): {len(string_fastq)} "
            "(Fun: fastq_to_dataframe)"
        )

        return string_fastq, df_fastq

    except Exception as e:
        logging.error(
            f"Failed to parse fastq file {fastq_file}: {e}. "
            "Ensure that the read2 fastq file is correctly formatted and accessible. "
            "You may also want to check if the file is corrupted or improperly indexed."
        )
        sys.exit(1)



#%%% Auxiliary functions to handle barcodes
#--------------------------------------------------------------------------------------------------------

def evaluate_barcode_quality(quality, threshold):

    """
    Evaluates the quality of a barcode based on its quality scores.

    Args:
        quality_scores (list): List of quality in ascii.
        threshold (int): Quality threshold.

    Returns:
        bool: True if the average quality meets or exceeds the threshold, False otherwise.
    """

    quality_scores = [ord(char) - 33 for char in quality]

    return len(quality_scores) > 0 and np.mean(np.array(quality_scores)) >= threshold

def compute_hamming_distance(my_bc, expected_bc):

    """
    Calculates the Hamming distance between two strings.

    Args:
        my_bc (str): Barcode to be evaluated.
        expected_bc (str): Expected barcode.

    Returns:
        int: Hamming distance, or None if the barcodes have different lengths.
    """
    if len(my_bc) != len(expected_bc):
        return None

    mybc_array = np.frombuffer(my_bc.encode(), dtype=np.uint8)
    expectedbc_array = np.frombuffer(expected_bc.encode(), dtype=np.uint8)

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

def correct_and_append(barcode, qbarcode, expected_data_dict, qthreshold, qname, bc_type):
    
    """
    Corrects the barcode if it matches the expected barcode set within a Hamming distance of 1.
    
    Args:
        barcode (str): The barcode to be corrected.
        qbarcode (list): Quality scores associated with the barcode.
        expected_data_dict (dict): Dictionary of expected barcodes with their well positions.
        qthreshold (int): Quality threshold to evaluate barcode quality.
        qname (str): Query name for logging purposes.
        bc_type (str): Type of barcode being corrected.

    Returns:
        str: The corrected barcode, or None if it can't be corrected.
    """

    if evaluate_barcode_quality(qbarcode, qthreshold):
        for expected_bc in expected_data_dict:
            hamming_dist = compute_hamming_distance(barcode, expected_bc)
            if hamming_dist == 1:
                return expected_bc  # Found a barcode with a Hamming distance of 1
            elif hamming_dist == 0:
                return barcode  # The barcode is already correct
    return None


#%%% Functions to obtain barcode coordinates
#--------------------------------------------------------------------------------------------------------

def run_bbmap(string_fastq, reference, num_threads, usr_memory_limit=0):

    """
    Runs bbmap.sh to align sequences and find barcodes.

    Args:
        fastq_as_df (pd.DataFrame): DataFrame containing BAM data.
        reference (str): Path to the reference genome.
        cores (int): Number of threads to use for parallel processing.
        usr_memory_limit (int): User-defined memory limit in MB (default is 0).

    Returns:
        bytes: Mapped BAM data in SAM format.
    """
    logging.info("Aligning sequence to find barcodes. This can take a while.")

    try:
        # Join the FASTQ entries into a single string and encode it
        fastq_data = '\n'.join(string_fastq).encode()

        # Calculate the memory limits
        available_memory = calculate_available_memory(usr_memory_limit)
        jvm_memory = calculate_jvm_memory(available_memory, safety_margin=0.7)

        # Run bbmap
        process = subprocess.Popen(
            [
                'bbmap.sh',
                f'ref={reference}',
                'in=stdin.fastq',
                'out=stdout.sam',
                f'-Xmx{jvm_memory}',
                'int=f',
                'k=8',
                'minid=0.1',
                'minratio=0.1',
                f'threads={num_threads}',
                'maxindel=100',
                'ambiguous=best',
                'local=false',
                'strictmaxindel=false'
             ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=False  # Important for binary input/output
        )

        bbmap_stdout, bbmap_stderr = process.communicate(input=fastq_data)

        logging.info("Sequence aligned to find barcodes")

        debug_logger.debug(
            f"bbmap output:\n\n"
            f"{bbmap_stderr.decode()}"
        )
            
        process.wait()

        return bbmap_stdout
    
    except subprocess.CalledProcessError as e:
        logging.error(
            f"Error occurred while running bbmap.sh: {e}. "
            "This may be due to issues with the reference genome, input FASTQ data, or system resources. "
            f"Check the stderr output for more details:\n{e.stderr.decode()}"
        )
        logging.error(
            "Ensure that the reference file exists, is correctly formatted, "
            "and that there is sufficient memory and CPU available."
        )
        sys.exit(1)

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

def parse_and_compute_coordinates(mapped_bam):

    """
    Parses the aligned BAM data to extract read names and CIGAR strings, 
    then computes the barcode coordinates based on CIGAR patterns.

    Args:
        mapped_bam (Data Frame): In-memory BAM data in SAM format.

    Returns:
        dict: Dictionary containing the computed barcode coordinates indexed by query names.
    """

    logging.info(f"Starting parsing barcodes to compute coordinates")

    records = []

    for read in mapped_bam[0]:

        if not read.strip():
            # Skip empty reads
            continue

        if read.startswith('@'):
            # Skip header
            continue

        fields = read.split('\t')
        cigar = fields[5]

        if cigar and re.match(r'^[0-9]+S.*[^0-9][7-9]I.*', cigar):
            records.append({
                'query_name': fields[0].split(' ')[0],  
                'cigarstring': cigar
                })

    map_df = pd.DataFrame(records)

    debug_logger.debug(
            f"Number of reads with barcode coordinates: {len(map_df)} "
            "(Fun: parse_and_compute_coordinates)"
        )

    # Get barcodes
    map_df = compute_initial_coordinates(map_df, 'BC1', r'([0-9]+S$)', 0)
    map_df = compute_initial_coordinates(map_df, 'BC2', r'([^0-9][7-9]I.*)', 0)
    map_df = compute_initial_coordinates(map_df, 'BC3', r'(?![0-9]).*$', -8)

    map_df = map_df[(map_df['BC1'] >= 17) & (map_df['BC2'] >= 9) & (map_df['BC3'] >= 1)]

    map_dict=map_df.set_index('query_name').to_dict('index')

    debug_logger.debug(
            f"Number of reads with corrected barcodes: {len(map_dict)} "
            "(Fun: parse_and_compute_coordinates)"
        )

    return map_dict
                        

#%%% Main function to process barcodes
#--------------------------------------------------------------------------------------------------------

def process_barcodes(clean_df, map_dict, barcode_files_dir, qthreshold, chunk_size=10000):
    """
    Processes barcode sequences from a DataFrame in chunks and applies barcode corrections, with a progress bar.

    Args:
        clean_df (pd.DataFrame): DataFrame containing sequence and quality information.
        map_dict (dict): Dictionary mapping query names to read mappings.
        barcode_files_dir (str): Directory containing expected barcode files.
        qthreshold (int): Quality threshold for barcode correction.
        chunk_size (int): Number of reads to process in each chunk.

    Returns:
        pd.DataFrame: Processed barcodes with corrected barcodes and UMI.
        pd.DataFrame: Aggregated statistics for barcode correction.
    """
    # Load barcode data into DataFrames
    bc1_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_1.csv"))
    bc2_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_2.csv"))
    bc3_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_3.csv"))

    # Create dictionaries for fast lookups
    bc1_dict = dict(zip(bc1_data['Barcode'], bc1_data['WellPosition']))
    bc2_dict = dict(zip(bc2_data['Barcode'], bc2_data['WellPosition']))
    bc3_dict = dict(zip(bc3_data['Barcode'], bc3_data['WellPosition']))

    barcode_ls = []
    barcode_stats_list = []

    for entry in clean_df.itertuples(index=False):
        if entry.query_name in map_dict:

            qname = entry.query_name
            rmap = map_dict[qname]
            seq = np.array(list(entry.seq))
            q = np.array(list(entry.qual))

            bc1 = ''.join(seq[rmap['BC1']:rmap['BC1'] + 10])
            bc2 = ''.join(seq[rmap['BC2']:rmap['BC2'] + 8])
            bc3 = ''.join(seq[rmap['BC3']:rmap['BC3'] + 8])

            qbc1 = ''.join(q[rmap['BC1']:rmap['BC1'] + 10])
            qbc2 = ''.join(q[rmap['BC2']:rmap['BC2'] + 8])
            qbc3 = ''.join(q[rmap['BC3']:rmap['BC3'] + 8])

            corrected_bc1 = correct_and_append(bc1, qbc1, bc1_dict, qthreshold, qname, 'Barcode1')
            corrected_bc2 = correct_and_append(bc2, qbc2, bc2_dict, qthreshold, qname, 'Barcode2')
            corrected_bc3 = correct_and_append(bc3, qbc3, bc3_dict, qthreshold, qname, 'Barcode3')

            umi_start = rmap['BC3'] - 11 if rmap['BC3'] > 10 else 0
            umi_end = rmap['BC3'] if rmap['BC3'] > 1 else 1
            umi = ''.join(seq[umi_start:umi_end])
            qumi = q[umi_start:umi_end]

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

            # Collect statistics for heatmaps
            for corrected_bc, bc_data in zip([corrected_bc1, corrected_bc2, corrected_bc3], [bc1_data, bc2_data, bc3_data]):
                if corrected_bc is not None:
                    well_position = bc_data.loc[bc_data['Barcode'] == corrected_bc, 'WellPosition'].values[0]
                    name = bc_data.loc[bc_data['Barcode'] == corrected_bc, 'Name'].values[0] if 'Name' in bc_data.columns else "Unknown"
                    barcode_stats_list.append({
                        'WellPosition': well_position,
                        'Name': name,
                        'Barcode': corrected_bc,
                        'Count': 1
                    })

    # Convert the list of barcodes and statistics to DataFrames
    barcode_df = pd.DataFrame(barcode_ls)

    debug_logger.debug(
            f"Number of reads with corrected barcodes: {len(barcode_df)} "
            "(Fun: process_barcodes)"
        )

    if barcode_stats_list:
        barcode_stats = pd.DataFrame(barcode_stats_list)
        barcode_stats = barcode_stats.groupby(['WellPosition', 'Name', 'Barcode'], as_index=False).sum()
    else:
        barcode_stats = pd.DataFrame(columns=['WellPosition', 'Name', 'Barcode', 'Count'])

    if barcode_stats.empty:
        debug_logger.debug(
                f"Barcode stats empty"
                "(Fun: process_barcodes)"
            )
    else:
        debug_logger.debug(
                f"Barcode stats are processed for {list(barcode_stats['Name'].unique())} "
                "(Fun: process_barcodes)"
            )
    
    return barcode_df, barcode_stats


#%%% Make graphs and stats
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

    # Check if data is empty
    if data.empty:
        logging.warning(f"No data available to create heatmap for {title}.")
        return

    heatmap_data = data.pivot_table(index="Row", columns="Column", values="Count", fill_value=0).infer_objects(copy=False)
    heatmap_data = heatmap_data.loc[sorted(heatmap_data.index, reverse=True)]

    indigo_cmap = LinearSegmentedColormap.from_list('white_to_indigo', ['white', '#002395'], N=1000)
    
    # Create the heatmap
    plt.figure(figsize=(12, 10))  
    ax = sns.heatmap(heatmap_data, annot=True, fmt=".0f", cmap=indigo_cmap, cbar=False, 
                     xticklabels=True, yticklabels=True, linewidths=0.5, linecolor='black')
    plt.title(title)
    plt.gca().invert_yaxis()
   
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.tick_params(left=False, top=False)
    
    plt.savefig(filename)
    plt.close()


#%%% Tag final bam
#--------------------------------------------------------------------------------------------------------

def tag_barcode(fastq_file, output_bam_file, corrected_barcode_df, chunk_size=1000):
    """
    Converts a FASTQ file to a BAM file and tags barcodes in chunks, with a progress bar.

    Args:
        fastq_file (str): Path to the input FASTQ file.
        output_bam_file (str): Path to the output BAM file.
        corrected_barcode_df (pd.DataFrame): DataFrame containing the corrected barcodes.
        chunk_size (int, optional): Number of reads to process in each chunk. Default is 10000.

    Returns:
        None
    """

    # Convert the DataFrame to a dictionary for quick lookup
    barcode_dict = corrected_barcode_df.set_index('query_name').to_dict('index')

    logging.info("Tagging barcodes in read1.")
    
    try:
        # Open the BAM file for writing
        with pysam.AlignmentFile(output_bam_file, "wb", header={'HD': {'VN': '1.0'}}) as bam_out:
            # Parse the FASTQ file
            with pysam.FastxFile(fastq_file) as fastq:
                total_reads = sum(1 for _ in pysam.FastxFile(fastq_file))  # Calculate the total number of entries
                with tqdm(total=total_reads, desc="Tagging barcodes") as pbar:
                    chunk = []
                    read_count = 0

                    for entry in fastq:
                        read_count += 1

                        # Create an unaligned BAM entry
                        read = pysam.AlignedSegment()
                        read.query_name = entry.name
                        read.query_sequence = entry.sequence
                        read.query_qualities = pysam.qualitystring_to_array(entry.quality)
                        read.flag = 4  # Unmapped flag

                        # Tag barcodes if available
                        if entry.name in barcode_dict:
                            barcode_info = barcode_dict[entry.name]
                            xd, xe, xf, xm = barcode_info['BC1'], barcode_info['BC2'], barcode_info['BC3'], barcode_info['UMI']
                            xc = xd + xe + xf

                            # Add tags to the read
                            read.set_tag('XD', xd, value_type='Z')
                            read.set_tag('XE', xe, value_type='Z')
                            read.set_tag('XF', xf, value_type='Z')
                            read.set_tag('XC', xc, value_type='Z')
                            read.set_tag('XM', xm, value_type='Z')

                        # Add the read to the chunk
                        chunk.append(read)

                        # Process the chunk if it reaches the specified size
                        if len(chunk) >= chunk_size:
                            for read in chunk:
                                bam_out.write(read)
                            chunk = []  # Reset the chunk

                        pbar.update(1)

                    # Write remaining reads in the chunk
                    if chunk:
                        for read in chunk:
                            bam_out.write(read)

        logging.info(f"Bam file with tagged barcodes in: {output_bam_file}")

    except Exception as e:
        logging.error(
            f"Failed to tag barcodes: {e}\n"
            "HELP - This might be due to an issue with read 1 fasta or memory limitations. "
            "Please check read1 fasta and system memory limitations.\n"
            "HELP - Use -v <verbose> for more detailed debugging"
        )
        sys.exit(1)


#%% PROCESS DATA
#########################################################################################################

def main():

    # Define the variables
    #----------------------------------------------------------------------------------------------------

    read1 = args.read1_fastq
    read2 = args.read2_fastq
    num_threads = args.threads
    barcode_files_dir = args.bc_dir
    qthreshold = args.qval
    output_bam = os.path.join(output_dir, f"{output_base}.bam")
    reference = os.path.join(barcode_files_dir, "invariable-linker-sequences.fasta")
    skip_tagging = args.skip_tagging
    usr_memory_limit= args.memory_limit * 1024
 
    
    # Open and parse fasta file for read 2
    #----------------------------------------------------------------------------------------------------
    
    logging.info("Starting barcode processing pipeline.")

    # Use ThreadPoolExecutor for opening and parsing the fastq file with the barcodes (read2)
    with ThreadPoolExecutor(max_workers=num_threads) as initial_executor:

        future_fastq = [initial_executor.submit(fastq_to_dataframe, read2)]

        clean_df = None
        string_fastq = None

        for future in as_completed(future_fastq):
            try:
                string_fastq, clean_df = future.result()
                logging.info("Parsed read2 to DataFrame")
                gc.collect()

            except Exception as e:
                logging.error(
                    f"Failed to parse read2 file: {e}. \n"
                    "HELP - Use -v <verbose> for more detailed debugging")
                sys.exit(1)   


    # Get coordinates
    #----------------------------------------------------------------------------------------------------

    if clean_df is not None:

        mapped_bam_data=run_bbmap(string_fastq, reference, num_threads, usr_memory_limit)
        logging.info("Barcode alignment to invariant sequence completed successfully.")

        # Adjusting memory limits to compute coordinates
        available_memory = calculate_available_memory(usr_memory_limit)
        num_processes = calculate_dynamic_process_count(available_memory, 
                                                        process_memory_overhead=500, 
                                                        safety_margin=0.7)
        mapped_bam_df=pd.DataFrame(mapped_bam_data.decode('utf-8').split('\n'))
        chunk_size = calculate_chunk_size(available_memory=available_memory, 
                                          memory_per_record=estimate_memory_per_record_df(mapped_bam_df))
        
        debug_logger.debug(
                f"Calculating memory use for parse_and_compute_coordinates:\n\n"
                f"Available memory: {int(available_memory // 1024)} GB\n"
                f"Number of parallel processes: {num_processes}\n"
                f"Chunk size: {chunk_size}\n"
            )
        
        total_records = len(mapped_bam_df)
        start_idx = 0

        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            futures_map_dict = []
            while start_idx < total_records:

                # Get the next chunk of data
                end_idx = min(start_idx + chunk_size, total_records)
                mapped_bam_chunk = mapped_bam_df.iloc[start_idx:end_idx]

                # Submit the chunk to ProcessPoolExecutor
                futures_map_dict.append(executor.submit(parse_and_compute_coordinates, mapped_bam_chunk))

                # Move to the next chunk
                start_idx = end_idx

            # Collect results
            for future_map_dict in as_completed(futures_map_dict):
                try:
                    map_dict = future_map_dict.result()
                    logging.info("Coordinates computed.")
                    
                except Exception as e:
                    logging.error(
                        f"Failed to align barcodes or parse coordinates: {e}.\n"
                        "HELP - Possible causes include issues with the reference file or improper barcode extraction. "
                        "Ensure all input files are correct.\n"
                        "HELP - Use -v <verbose> for more detailed debugging"
                    )
                    sys.exit(1)

                finally:
                    gc.collect()  
    
    # Correct barcodes
    #----------------------------------------------------------------------------------------------------

    original_reads = len(clean_df) if clean_df is not None else 0
    tagged_reads = 0
    efficiency = 0.0
    num_cells = 0

    if map_dict and not clean_df.empty:

        # Adjusting memory limits to compute coordinates
        available_memory = calculate_available_memory(usr_memory_limit)
        num_processes = calculate_dynamic_process_count(available_memory, 
                                                        process_memory_overhead=500, 
                                                        safety_margin=0.7)
        chunk_size = calculate_chunk_size(available_memory=available_memory, 
                                          memory_per_record=estimate_memory_per_record_df(clean_df))
        
        debug_logger.debug(
                f"Calculating memory use for process_barcodes:\n\n"
                f"Available memory: {int(available_memory // 1024)} GB\n"
                f"Number of parallel processes: {num_processes}\n"
                f"Chunk size: {chunk_size}\n"
            )
        
        total_records = len(clean_df)
        start_idx = 0

        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            futures_processedbc = []
            while start_idx < total_records:

                # Get the next chunk of data
                end_idx = min(start_idx + chunk_size, total_records)
                clean_df_chunk = clean_df.iloc[start_idx:end_idx]

                # Submit the chunk to ProcessPoolExecutor
                futures_processedbc.append(executor.submit(process_barcodes, clean_df_chunk, map_dict, barcode_files_dir=barcode_files_dir, qthreshold=qthreshold))

                # Move to the next chunk
                start_idx = end_idx

            # Collect results
            for future_processedbc in as_completed(futures_processedbc):
                try:
                    barcode_df, stats_df = future_processedbc.result()
                    logging.info("Barcodes corrected")
                    
                except Exception as e:
                    logging.error(
                        f"Failed to process barcodes: {e}\n"
                        "HELP - This might be due to an issue with the input data, such as not finding the barcodes. "
                        "Please check the input files and the barcode directory\n"
                        "HELP - Use -v <verbose> for more detailed debugging"
                    )
                    sys.exit(1)

                finally:
                    gc.collect()  


    # Calculate stats and crate heatmaps
    #----------------------------------------------------------------------------------------------------

    if not barcode_df.empty:

        tagged_reads = len(barcode_df)
        efficiency = (tagged_reads / original_reads) * 100 if original_reads > 0 else 0.0

        cell_counts = barcode_df.groupby(['BC1', 'BC2', 'BC3']).size().reset_index(name='counts')
        valid_cells = cell_counts[cell_counts['counts'] > 10]
        num_cells = len(valid_cells)

        parsed_df = parse_well_position(stats_df)
        for barcode in parsed_df['Name'].unique():
            barcode_data = parsed_df[parsed_df['Name'] == barcode]
            heatmap_filename = os.path.join(output_dir, f"{output_base}_heatmap_{barcode}.png")
            logging.info(f"Creating heatmap for {barcode} in: {heatmap_filename}")
            create_heatmap(barcode_data, f"Heatmap for {barcode}", heatmap_filename)
                    
            gc.collect()

    else:
        logging.error(
            f"Processed barcode DataFrame is empty.\n"
            "HELP - This might be due to an issue with the input data, such as not finding the barcodes. "
            "Please check the input files and the barcode directory\n"
            "HELP - Use -v <verbose> for more detailed debugging"
            )
        sys.exit(1)


    # Tag barcodes in a bam file (if asked)
    #----------------------------------------------------------------------------------------------------

    if not skip_tagging:

        # Adjusting memory limits to tag_barcode
        available_memory = calculate_available_memory(usr_memory_limit)
        chunk_size = calculate_chunk_size(available_memory=available_memory, 
                                          memory_per_record=estimate_memory_per_record_fq(read1),
                                          threads=num_threads)
        debug_logger.debug(
                f"Calculating memory use for tag_barcode:\n\n"
                f"Available memory: {int(available_memory // 1024)} GB\n"
                f"Number used threads: {num_threads}\n"
                f"Chunk size: {chunk_size}\n"
            )
        
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            executor.submit(tag_barcode, 
                            fastq_file=read1,
                            output_bam_file=output_bam,
                            corrected_barcode_df=barcode_df,
                            chunk_size=chunk_size)

    else:
        logging.info("Skipping the tagging step as requested.")

    
    # Print stats_df to stats_logger (if --stats is inputed)
    #----------------------------------------------------------------------------------------------------

    if stats_logger:
        stats_df = stats_df.sort_values(by=['Name', 'WellPosition'])
        stats_df_dropped = stats_df.drop(columns=['Column', 'Row'])
        stats_logger.info(stats_df_dropped.to_string(index=False))

    # Finileze script
    #----------------------------------------------------------------------------------------------------

    logging.info("Barcode processing pipeline completed.")
    logging.info("Have fun analyzing! (>ᴗ•)❀\n")

    # Print final statistics in the log file
    logging.info(
        f"Final info:\n"
        f"\t\t\tNumber of original reads: {original_reads}\n"
        f"\t\t\tNumber of taggable reads: {tagged_reads}\n"
        f"\t\t\tEfficiency: {efficiency:.2f}%\n"
        f"\t\t\tNumber of unique cells: {num_cells}"
    )

    # Ensure all logging is outputted
    for handler in logging.root.handlers[:]:
        handler.flush()
        handler.close()


if __name__ == "__main__":
    main()