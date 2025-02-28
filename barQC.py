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

# Standard Library Imports
import argparse
import os
import warnings
import logging
import subprocess
import sys
import gc
import tempfile
import re
import psutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import pysam
from pysam import FastxFile
from concurrent.futures import ThreadPoolExecutor, as_completed

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
optional.add_argument("--temp_dir", 
                      help="Optional directory to store temporary files. Defaults to system temp directory.",
                      default=None)

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
debug_log_path = os.path.join(output_dir, f"{output_base}_debug.log")
debug_logger = logging.getLogger('debug_logger')

if args.verbose:
    debug_logger.setLevel(logging.DEBUG)
    debug_handler = logging.FileHandler(debug_log_path, mode='w')
    debug_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    debug_logger.addHandler(debug_handler)
    debug_logger.propagate = False 
else:
    debug_logger.setLevel(logging.WARNING)
    debug_logger.addHandler(logging.NullHandler())


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

 # Validate if temp_dir exists
if args.temp_dir and not os.path.isdir(args.temp_dir):
    print(f"\n\nERROR: Invalid temp_dir. Directory does not exist: {args.temp_dir}\n")
    logging.error(f"Invalid temp_dir. Directory does not exist: {args.temp_dir}\n")
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
        if args.verbose:
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
    """
    Estimates memory usage per record for a FASTQ file.

    Args:
        file_path (str): Path to the FASTQ file.

    Returns:
        tuple: Number of records (int) and memory usage per record in MB (float).
    """
    file_size_bytes = os.path.getsize(file_path)

    with open(file_path, 'r') as f:
        num_lines = sum(1 for _ in f)
    num_records = num_lines // 4

    if num_records == 0:
        logging.error(f"Unable to process FASTA file\n"
            "HELP - Use -v <verbose> for more detailed debugging")
        debug_logger.debug("Unable to estimate memory usage per record. File might be empty or improperly formatted.")
        sys.exit(1)

    memory_per_record_bytes = file_size_bytes / num_records
    memory_per_record_mb = memory_per_record_bytes / (1024 * 1024)

    return num_records,memory_per_record_mb

def calculate_available_memory(usr_memory_limit):
    """
    Calculates available memory based on user-defined limits and system capacity.

    Args:
        usr_memory_limit (int): User-specified memory limit in MB.

    Returns:
        int: Available memory in MB.
    """

    if usr_memory_limit > 0:
        system_memory = psutil.virtual_memory().available // (1024 * 1024)
        available_memory = min(usr_memory_limit, system_memory)
    else:
        available_memory = psutil.virtual_memory().available // (1024 * 1024)

    return available_memory

def calculate_jvm_memory(available_memory, safety_margin):
    """
    Calculates the JVM memory allocation based on available memory and a safety margin.

    Args:
        available_memory (int): Available memory in MB.
        safety_margin (float): Fraction of available memory to allocate.

    Returns:
        int: JVM memory allocation in GB.
    """


    jvm_memory_mb = int(available_memory * safety_margin)
    
    jvm_memory_gb = jvm_memory_mb // 1024
    return jvm_memory_gb

def calculate_chunk_size(available_memory, memory_per_record, safety_margin, threads):
    """
    Calculates the optimal chunk size for processing records.

    Args:
        available_memory (int): Available memory in MB.
        memory_per_record (float): Memory required per record in MB.
        safety_margin (float): Fraction of memory to safely use.
        threads (int): Number of processing threads.

    Returns:
        int: Calculated chunk size
    """

    usable_memory = (available_memory * safety_margin) / threads
    chunk_size = int(usable_memory / memory_per_record)
    return chunk_size


#%%% Auxiliary functions to handle barcodes
#--------------------------------------------------------------------------------------------------------

def evaluate_barcode_quality(quality, threshold):
    """
    Evaluates the quality of a barcode based on its quality scores.

    Args:
        quality (str): String of ASCII-encoded quality scores.
        threshold (float): Minimum average quality score for acceptance.

    Returns:
        bool: True if the barcode meets the quality threshold, False otherwise.
    """

    quality_scores = [ord(char) - 33 for char in quality]

    return len(quality_scores) > 0 and np.mean(np.array(quality_scores)) >= threshold

def compute_hamming_distance(my_bc, expected_bc):
    """
    Computes the Hamming distance between two barcodes.

    Args:
        my_bc (str): Observed barcode.
        expected_bc (str): Expected barcode.

    Returns:
        int or None: Hamming distance if barcodes are of equal length, otherwise None.
    """


    if len(my_bc) != len(expected_bc):
        return None

    mybc_array = np.frombuffer(my_bc.encode(), dtype=np.uint8)
    expectedbc_array = np.frombuffer(expected_bc.encode(), dtype=np.uint8)

    return np.sum(mybc_array != expectedbc_array)

def correct_and_append(barcode, qbarcode, expected_data_dict, qthreshold):
    """
    Computes the Hamming distance between two barcodes.

    Args:
        my_bc (str): Observed barcode.
        expected_bc (str): Expected barcode.

    Returns:
        int or None: Hamming distance if barcodes are of equal length, otherwise None.
    """

    if not evaluate_barcode_quality(qbarcode, qthreshold):
        return None

    best_match = None
    min_distance = 3  # Only accept corrections within Hamming distance of 2

    for expected_bc in expected_data_dict:
        if barcode == expected_bc:
            return expected_bc

        hamming_dist = compute_hamming_distance(barcode, expected_bc)
        if hamming_dist is not None and hamming_dist < min_distance:
            best_match = expected_bc
            min_distance = hamming_dist
            if hamming_dist == 1:
                break

    return best_match

def calculate_barcode_length(df):
    """
    Calculates the length of barcodes in a DataFrame, ensuring uniformity.

    Args:
        df (pandas.DataFrame): DataFrame containing barcodes in the third column.

    Returns:
        int: Length of barcodes if all are uniform.

    Raises:
        Logs an error if barcodes have varying lengths.
    """

    lengths = df.iloc[:, 2].apply(len)
    if lengths.nunique() == 1:
        return lengths.iloc[0]
    else:
        raise ValueError("Expected barcodes have different sizes")


#%%% Auxiliary functions to make graphs and stats
#--------------------------------------------------------------------------------------------------------

def parse_well_position(df):
    """
    Parses well positions into row and column components.

    Args:
        df (pandas.DataFrame): DataFrame containing a 'WellPosition' column.

    Returns:
        pandas.DataFrame: Updated DataFrame with 'Row' and 'Column' columns added.
    """

    df['Row'] = df['WellPosition'].str[0]
    df['Column'] = df['WellPosition'].str[1:].astype(int)
    return df

def create_heatmap(data, title, filename):
    """
    Creates and saves a heatmap from the provided data.

    Args:
        data (pandas.DataFrame): DataFrame with 'Row', 'Column', and 'Count' columns.
        title (str): Title of the heatmap.
        filename (str): File path to save the heatmap image.

    Returns:
        None: Saves the heatmap to a file and exits the function.
    """

    # Check if data is empty
    if data.empty:
        logging.warning(f"No data available to create heatmap for {title}.")
        return

    heatmap_data = data.pivot_table(index="Row", columns="Column", values="Count", fill_value=0).infer_objects(copy=False)
    heatmap_data = heatmap_data.loc[sorted(heatmap_data.index, reverse=True)]

    indigo_cmap = LinearSegmentedColormap.from_list('white_to_indigo', ['white', '#1e568b'], N=1000)
    
    # Make the ploting values in log to have a better scale
    heatmap_data_log = np.log10(heatmap_data + 1)

    # Create the heatmap
    # Define figure size and initialize heatmap
    plt.figure(figsize=(12, 8))

    # Create heatmap
    ax = sns.heatmap(
        heatmap_data_log,
        annot=heatmap_data,
        fmt=".0f",
        cmap=indigo_cmap,
        cbar=False,
        xticklabels=True,
        yticklabels=True,
        linewidths=0.5,
        linecolor='black'
    )

    # Add title with custom styling
    plt.title(
        title,
        fontsize=16,
        fontweight='bold',
        pad=35
    )

    # Invert y-axis for better visual alignment
    plt.gca().invert_yaxis()

    # Remove default axis labels
    plt.xlabel('')
    plt.ylabel('')

    # Adjust tick label orientation
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)

    # Move x-axis labels to the top
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

    # Remove unnecessary tick marks
    ax.tick_params(left=False, top=False)

    # Adjust layout to prevent clipping
    plt.tight_layout()

    # Save the figure to file & close plot
    plt.savefig(filename, bbox_inches='tight')
    plt.close()


#%%% Functions handel fasta files
#--------------------------------------------------------------------------------------------------------

def split_pairfasta_to_temp(read1, read2, temp_dir, memory_limit, threads, safety_margin=0.85, max_chunk_size=100000):
    """
    Splits paired FASTA files into smaller chunks for parallel processing.

    Args:
        read1 (str): Path to the first paired FASTA file.
        read2 (str): Path to the second paired FASTA file.
        temp_dir (str): Directory to store temporary chunk files.
        memory_limit (int): Memory limit in MB.
        threads (int): Number of processing threads.
        safety_margin (float, optional): Fraction of memory to safely allocate. Default is 0.85.
        max_chunk_size (int, optional): Maximum records per chunk. Default is 100,000.

    Returns:
        tuple: List of created temporary file pairs and the total number of records.

    Raises:
        Logs errors for issues during processing and performs cleanup before exiting.
    """

    temp_files = []

    try:
        # Estimate memory usage per record
        num_records,memory_per_record_mb = estimate_memory_per_record_fq(read2)
        
        # Calculate available memory and dynamic chunk size
        available_memory = calculate_available_memory(memory_limit)
        calculated_chunk_size = calculate_chunk_size(available_memory, memory_per_record_mb, safety_margin, threads)

        # Force a maximum chunk size for better parallelization
        chunk_size = min(calculated_chunk_size, max_chunk_size)
        chunk_size = (chunk_size // 4) * 4

        debug_logger.debug(f"You are working with {num_records} reads")
        debug_logger.debug(f"Estimated memory per record: {memory_per_record_mb:.2f} MB")
        debug_logger.debug(f"Available memory: {available_memory} MB")
        debug_logger.debug(f"Temporal fasta size: {chunk_size} records")

        with FastxFile(read1) as fq1, FastxFile(read2) as fq2:
            while True:
                paired_records = []

                for _ in range(chunk_size):
                    record1 = next(fq1, None)
                    record2 = next(fq2, None)

                    if record1 is None or record2 is None:
                        break

                    paired_records.append((record1, record2))

                if not paired_records:
                    break
                
                # Create temporary files for the current chunk
                temp_file1 = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fastq', dir=temp_dir)
                temp_file2 = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fastq', dir=temp_dir)
                
                # Write records to temporary files
                for record1, record2 in paired_records:
                    if record1 and record2:  # Ensure both reads exist
                        temp_file1.write(f"@{record1.name}\n{record1.sequence}\n+\n{record1.quality}\n")
                        temp_file2.write(f"@{record2.name}\n{record2.sequence}\n+\n{record2.quality}\n")
                
                temp_file1.close()
                temp_file2.close()
                temp_files.append((temp_file1.name, temp_file2.name))

        logging.info(
            f"Created {len(temp_files)} temporary files for parallel processing"
        )
        return temp_files, num_records

    except Exception as e:
        logging.error(f"Error while splitting FASTA file\n"
            "HELP - Use -v <verbose> for more detailed debugging")
        debug_logger.debug(f"Error while splitting FASTA file: {e}")

        # Clean up created files if an error occurs
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        logging.info(f"The program encountered a critical issue and is exiting.")
        sys.exit(1)

def split_fasta2_to_temp(read2, temp_dir, memory_limit, threads, safety_margin=0.85, max_chunk_size=100000):
    """
    Splits a FASTA file into temporary files for parallel processing, dynamically adjusting chunk size based on memory constraints.

    Args:
        read2 (str): Path to the input FASTA file.
        temp_dir (str): Directory to store temporary chunk files.
        memory_limit (int): User-specified memory limit in MB.
        threads (int): Number of threads for parallel processing.
        safety_margin (float, optional): Fraction of memory to safely allocate. Default is 0.85.
        max_chunk_size (int, optional): Maximum number of records per chunk. Default is 100,000.

    Returns:
        tuple: A list of paths to created temporary files and the total number of records.

    Raises:
        Logs an error if processing fails, cleans up temporary files, and exits the program.
    """

    temp_files = []

    try:
               # Estimate memory usage per record
        num_records,memory_per_record_mb = estimate_memory_per_record_fq(read2)
        
        # Calculate available memory and dynamic chunk size
        available_memory = calculate_available_memory(memory_limit)
        calculated_chunk_size = calculate_chunk_size(available_memory, memory_per_record_mb, safety_margin, threads)

        # Force a maximum chunk size for better parallelization
        chunk_size = min(calculated_chunk_size, max_chunk_size)
        chunk_size = (chunk_size // 4) * 4

        debug_logger.debug(f"You are working with {num_records} reads")
        debug_logger.debug(f"Estimated memory per record: {memory_per_record_mb:.2f} MB")
        debug_logger.debug(f"Available memory: {available_memory} MB")
        debug_logger.debug(f"Temporal fasta size: {chunk_size} records")

        with FastxFile(read2) as fq2:
            while True:
                records = []

                for _ in range(chunk_size):
                    record = next(fq2, None)

                    if record is None:
                        break

                    records.append(record)

                if not records:
                    break
                
                # Create a temporary file for the current chunk
                temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fastq', dir=temp_dir)
                
                # Write records to the temporary file
                for record in records:
                    temp_file.write(f"@{record.name}\n{record.sequence}\n+\n{record.quality}\n")
                
                temp_file.close()
                temp_files.append(temp_file.name)

        logging.info(
            f"Created {len(temp_files)} temporary files for parallel processing."
        )
        return temp_files, num_records

    except Exception as e:
        logging.error(f"Error while splitting FASTA file\n"
            "HELP - Use -v <verbose> for more detailed debugging")
        debug_logger.debug(f"Error while splitting FASTA file: {e}")

        # Clean up created files if an error occurs
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        logging.info(f"The program encountered a critical issue and is exiting.")
        sys.exit(1)

def truncate_fastq(fastq, output_fastq, length=100):
    """
    Truncates each sequence and quality in a FASTQ file to a specified length.

    Args:
        fastq (str): Path to the input FASTQ file.
        output_fastq (str): Path to the output truncated FASTQ file.
        length (int, optional): Number of bases to keep for each sequence. Default is 100.

    Returns:
        None: Writes the truncated FASTQ entries to the output file.
    """

    with open(fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        while True:
            # Read four lines for each FASTQ entry
            header = infile.readline()
            if not header:  # End of file
                break
            sequence = infile.readline().strip()
            plus = infile.readline()
            quality = infile.readline().strip()
            
            # Truncate sequence and quality
            truncated_sequence = sequence[:length]
            truncated_quality = quality[:length]
            
            # Write truncated entry
            outfile.write(header)
            outfile.write(truncated_sequence + '\n')
            outfile.write(plus)
            outfile.write(truncated_quality + '\n')

#%%% Functions to obtain barcode coordinates & correct
#--------------------------------------------------------------------------------------------------------

def run_bbmap(fastq, reference, num_threads, usr_memory_limit):


    try:

        # Create a temporary file for truncated FASTQ
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fastq") as temp_fastq:
            truncated_fastq = temp_fastq.name
        
        # Truncate the FASTQ file to the first 100 bases
        truncate_fastq(fastq, truncated_fastq, length=100)

        # Calculate the memory limits
        available_memory = calculate_available_memory(usr_memory_limit)
        jvm_memory = calculate_jvm_memory(available_memory, safety_margin=0.85)

        # Run bbmap
        process = subprocess.Popen(
            [
                'bbmap.sh',
                f'ref={reference}',
                f'in={truncated_fastq}',
                'out=stdout.sam',
                f'-Xmx{max(jvm_memory - 2, 4)}g',
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

        bbmap_stdout, bbmap_stderr = process.communicate()


        decoded_stderr = bbmap_stderr.decode() 
        debug_logger.debug(
            f"bbmap run:\n"
            f"{'\n'.join(decoded_stderr.splitlines()[:3])}\n"
            f"{decoded_stderr.splitlines()[-1]}"
        )
        
        # Clean up the temporary file
        os.remove(truncated_fastq)

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

    records = []

    for read in mapped_bam[0]:
        if not read.strip():
            continue

        if read.startswith('@'):
            continue

        fields = read.split('\t')
        if len(fields) <= 5:
            continue
        
        cigar = fields[5]
        if not cigar or not re.match(r'^[0-9]+S.*[^0-9][7-9]I.*', cigar):
            continue

        print(records)
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

    map_dict=map_df.set_index('query_name').to_dict('index')

    return map_dict
         
def process_barcodes(read2, map_dict, expectedbc_df_list, expectedbc_len_list, qthreshold):

    # Load barcode data into dictionaries
    dict_list = []
    for _, df in enumerate(expectedbc_df_list):
        bc_dict = dict(zip(df['Barcode'], df['WellPosition']))
        dict_list.append(bc_dict)
    
    bc1_dict = dict_list[0]
    bc2_dict = dict_list[1]
    bc3_dict = dict_list[2]

    barcode_ls = []

    try:
        with pysam.FastxFile(read2) as fastq:
            for entry in fastq:
                if entry.name in map_dict:

                    qname = entry.name
                    rmap = map_dict[qname]
                    seq = np.array(list(entry.sequence))
                    q = np.array(list(entry.quality))

                    bc1 = ''.join(seq[rmap['BC1']:rmap['BC1'] + expectedbc_len_list[0]])
                    bc2 = ''.join(seq[rmap['BC2']:rmap['BC2'] + expectedbc_len_list[1]])
                    bc3 = ''.join(seq[rmap['BC3']:rmap['BC3'] + expectedbc_len_list[2]])

                    qbc1 = ''.join(q[rmap['BC1']:rmap['BC1'] + expectedbc_len_list[0]])
                    qbc2 = ''.join(q[rmap['BC2']:rmap['BC2'] + expectedbc_len_list[1]])
                    qbc3 = ''.join(q[rmap['BC3']:rmap['BC3'] + expectedbc_len_list[2]])

                    corrected_bc1 = correct_and_append(bc1, qbc1, bc1_dict, qthreshold)
                    corrected_bc2 = correct_and_append(bc2, qbc2, bc2_dict, qthreshold)
                    corrected_bc3 = correct_and_append(bc3, qbc3, bc3_dict, qthreshold)

                    umi_start = rmap['BC3'] - 11 if rmap['BC3'] > 10 else 0
                    umi_end = rmap['BC3'] if rmap['BC3'] > 1 else 1
                    umi = ''.join(seq[umi_start:umi_end])
                    qumi = q[umi_start:umi_end]

                    if len(umi) != 10 or not evaluate_barcode_quality(qumi, qthreshold):
                        continue

                    if None not in [corrected_bc1, corrected_bc2, corrected_bc3] and umi is not None:
                        barcode_ls.append({
                            'query_name': qname,
                            'BC1': corrected_bc1,
                            'BC2': corrected_bc2,
                            'BC3': corrected_bc3,
                            'UMI': umi
                        })

    except Exception as e:

        logging.error(
            f"Failed to parse fastq file {read2}\n"
            "HELP - Use -v <verbose> for more detailed debugging")
        logging.info(
            "Ensure that the read2 fastq file is correctly formatted and accessible. "
            "You may also want to check if the file is corrupted"
        )
        debug_logger.debug(
            f"Failed to parse fastq file {read2}: {e}")
        
        sys.exit(1)

    # Convert the list of barcodes and statistics to DataFrames
    barcode_df = pd.DataFrame(barcode_ls)
    return barcode_df

#%% Functions to obteing stats
#--------------------------------------------------------------------------------------------------------

def calculate_stats(barcode_df, expectedbc_df_list):

    try:
        # Barcode stats
        #--------------

        # Validate barcode_df
        if barcode_df.empty:
            logging.error("No statistics or heatmaps generated.\n"
                    "HELP - Use -v <verbose> for more detailed debugging")
            sys.exit(1)

        # Create the final DataFrame to hold results
        barcode_stats_list = []

        # Loop through the expected barcode DataFrames
        for i, bc_df in enumerate(expectedbc_df_list):
            col_name = f"BC{i+1}"
            for _, row in bc_df.iterrows():
                barcode = row['Barcode']
                count = barcode_df[col_name].value_counts().get(barcode, 0)
                barcode_stats_list.append({
                    'WellPosition': row['WellPosition'],
                    'Name': row['Name'],
                    'Barcode': barcode,
                    'Count': count
                })

        # Convert results to a DataFrame
        barcode_stats = pd.DataFrame(barcode_stats_list)

        # UMI stats
        #----------
        umi_stats = pd.DataFrame()
        umi_stats['Cell'] = barcode_df['BC1'] + barcode_df['BC2'] + barcode_df['BC3']
        umi_stats['MolecID'] = barcode_df['BC1'] + barcode_df['BC2'] + barcode_df['BC3'] + barcode_df['UMI']
        
        umi_stats['UMI_Counts'] = umi_stats.groupby("Cell")["MolecID"].transform('count')
        umi_stats = umi_stats.drop_duplicates()

        return barcode_stats, umi_stats
    
    except Exception as e:
        logging.error(f"Error calculating statistics\n"
                    "HELP - Use -v <verbose> for more detailed debugging")
        debug_logger.debug(f"Error calculating statistics: {e}")
    
def create_heatmaps(barcode_stats, output_dir, output_base):
    
    try:    
        barcode_dfs = {name: group.reset_index(drop=True) for name, group in barcode_stats.groupby("Name")}

        for bcname, data in barcode_dfs.items():
            parsed_df = parse_well_position(data)
            heatmap_filename = os.path.join(output_dir, f"{output_base}_heatmap_{bcname}.png")
            create_heatmap(parsed_df, 
                           title=f"Heatmap for {bcname}", 
                           filename=heatmap_filename)
    except Exception as e:
        logging.error(f"Error creating heatmaps\n"
                    "HELP - Use -v <verbose> for more detailed debugging")
        debug_logger.debug(f"Error creating heatmaps: {e}")

def create_umi_plots(umi_stats, output_dir, output_base):
    
    # Plot UMI duplicaction
    plt.figure(figsize=(12, 6))
    umi_stats['UMI_Counts'].plot(kind='hist', 
                                 bins=range(1, umi_stats['UMI_Counts'].max() + 2), 
                             color='#2E6BA1', edgecolor='#0A1C2E', 
                             alpha=0.9, align='left')
    plt.title("Distribution of UMI Counts", fontsize=16, fontweight='bold')
    plt.xlabel("UMI Counts", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{output_base}_UMI_counts_distribution.png"))
    plt.close()

    
    # Plot mean UMI duplication per cell
    cell_mean_counts = umi_stats.groupby('Cell')['UMI_Counts'].mean()
    plt.figure(figsize=(12, 6))
    cell_mean_counts.plot(kind='hist', 
                      bins=np.linspace(cell_mean_counts.min(), cell_mean_counts.max(), 20),  # 20 bins
                      color='#2E6BA1', edgecolor='#0A1C2E', 
                      alpha=0.9)
    plt.title("Distribution of Mean UMI Counts per Unique Barcode Combination", fontsize=16, fontweight='bold')
    plt.xlabel("Mean UMI Counts", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{output_base}_UMI_percell_distribution.png"))
    plt.close()
    

#%%% Functions to create and tag the bam file
#--------------------------------------------------------------------------------------------------------

def tag_barcode(fastq_file, output_bam_file, corrected_barcode_df, header):

    # Convert the DataFrame to a dictionary for quick lookup
    barcode_dict = corrected_barcode_df.set_index('query_name').to_dict('index')
    
    try:
        # Open the BAM file for writing
        with pysam.AlignmentFile(output_bam_file, "wb", header=header) as bam_out:
            
            # Parse the FASTQ file
            with pysam.FastxFile(fastq_file) as fastq:
                for entry in fastq:
                    # Process only reads in barcode_dict
                    if entry.name in barcode_dict:
                        # Create an unaligned BAM entry
                        read = pysam.AlignedSegment()
                        read.query_name = entry.name
                        read.query_sequence = entry.sequence
                        read.query_qualities = pysam.qualitystring_to_array(entry.quality)
                        read.flag = 4  # Unmapped flag

                        # Tag barcodes if available
                        barcode_info = barcode_dict[entry.name]
                        xd, xe, xf, xm = barcode_info['BC1'], barcode_info['BC2'], barcode_info['BC3'], barcode_info['UMI']
                        xc = xd + xe + xf

                        # Add tags to the read
                        read.set_tag('XD', xd, value_type='Z')
                        read.set_tag('XE', xe, value_type='Z')
                        read.set_tag('XF', xf, value_type='Z')
                        read.set_tag('XC', xc, value_type='Z')
                        read.set_tag('XM', xm, value_type='Z')

                        # Write the read to the BAM file
                        bam_out.write(read)

    except Exception as e:
        logging.error(
            f"Failed to tag barcodes: {e}\n"
            "HELP - This might be due to an issue with read1 FASTA or memory limitations. "
            "Please check read1 FASTA and system memory limitations.\n"
            "HELP - Use -v <verbose> for more detailed debugging"
        )
        sys.exit(1)


#%% Barcode correction wrappers
#--------------------------------------------------------------------------------------------------------

def barcode_correction_pipeline(read2, barcode_files_dir, reference, qthreshold, num_threads, usr_memory_limit):

    # Get coordinates
    #----------------
    mapped_bam_data = run_bbmap(read2, reference, num_threads, usr_memory_limit)
    mapped_bam_df = pd.DataFrame(mapped_bam_data.decode('utf-8').split('\n'))

    try:
        map_dict = parse_and_compute_coordinates(mapped_bam_df)
    
    except Exception as e:
        logging.error(
            f"Failed to align barcodes or parse coordinates: {e}.\n"
            "HELP - Possible causes include issues with the reference file or improper barcode extraction. "
            "Ensure all input files are correct.\n"
            "HELP - Use -v <verbose> for more detailed debugging")
        sys.exit(1)

    finally:
        gc.collect()
    

    # Correct barcodes
    #-----------------
    try:
        barcode_df = process_barcodes(read2, map_dict, expectedbc_df_list, expectedbc_len_list,qthreshold)
        # print(barcode_df)
    
    except Exception as e:
        logging.error(
            f"Failed to process barcodes\n"
            "HELP - Use -v <verbose> for more detailed debugging")
        logging.info(
            "HELP - This might be due to an issue with the input data, such as not finding the barcodes. "
            "Please check the input files and the barcode directory")

        debug_logger.debug(
            f"Failed to process barcodes: {e}")
        sys.exit(1)
    
    finally:
        gc.collect()

    return barcode_df

def concatenate_bam_files(temp_bam_files, output_bam):
    """
    Concatenates multiple BAM files into a single BAM file using samtools.

    Args:
        temp_bam_files (list): List of paths to temporary BAM files.
        output_bam (str): Path to the final merged BAM file.

    Returns:
        str: Path to the concatenated BAM file.
    """
    try:
        # Build samtools merge command
        merge_command = ["samtools", "merge", "-f", output_bam] + temp_bam_files
        subprocess.run(merge_command, check=True)
        return output_bam
    except Exception as e:
        logging.error(f"Error while merging BAM files: {e}")
        raise
    finally:
        gc.collect()

#%% PROCESS DATA
#########################################################################################################


#%% Define the variables
#----------------------------------------------------------------------------------------------------
read1 = args.read1_fastq
read2 = args.read2_fastq
num_threads = args.threads
barcode_files_dir = args.bc_dir
temp_dir = args.temp_dir
qthreshold = args.qval
output_bam = os.path.join(output_dir, f"{output_base}.bam")
reference = os.path.join(barcode_files_dir, "invariable-linker-sequences.fasta")
skip_tagging = args.skip_tagging
usr_memory_limit= args.memory_limit * 1024
bam_header = {'HD': {'VN': '1.0'},
              'PG': [{'ID': 'barQC',
                      'VN': '1.0',
                      'CL': f"barQC.py -f1 {read1} -f2 {read2} -o {output_bam} "
                            f"-t {num_threads} -d {barcode_files_dir} -q {qthreshold} "
                            f"--memory_limit {usr_memory_limit // 1024} "
                            f"{'--skip_tagging' if skip_tagging else ''}"}]}

#%% Load barcode files
#----------------------------------------------------------------------------------------------------
bc1_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_1.csv"))
bc2_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_2.csv"))
bc3_data = pd.read_csv(os.path.join(barcode_files_dir, "expected_barcodes_3.csv"))
expectedbc_df_list = [bc1_data, bc2_data, bc3_data]

try:
    bc1_length = calculate_barcode_length(bc1_data)
    debug_logger.debug(f"Barcode 1 has {bc1_length} bases")

except ValueError as e:
    logging.error(f"Error in Barcode 1 expected barcodes: {e}\n"
    "HELP - Check thet all expected barcodes in Barcode 1 have the same size and are correctly formated\n"
    "Expected tab separated table with headers: WellPosition Name Barcode")
    sys.exit(1)

try:
    bc2_length = calculate_barcode_length(bc2_data)
    debug_logger.debug(f"Barcode 2 has {bc2_length} bases")
    
except ValueError as e:
    logging.error(f"Error in Barcode 2 expected barcodes: {e}\n"
    "HELP - Check thet all expected barcodes in Barcode 2 have the same size and are correctly formated\n"
    "Expected tab separated table with headers: WellPosition Name Barcode")
    sys.exit(1)

try:
    bc3_length = calculate_barcode_length(bc3_data)
    debug_logger.debug(f"Barcode 3 has {bc3_length} bases")
    
except ValueError as e:
    logging.error(f"Error in Barcode 3 expected barcodes: {e}\n"
    "HELP - Check thet all expected barcodes in Barcode 3 have the same size and are correctly formated\n"
    "Expected tab separated table with headers: WellPosition Name Barcode")
    sys.exit(1)

expectedbc_len_list = [bc1_length, bc2_length, bc3_length]

logging.info(f"Starting to run BarQC!")

#%% Split the input in smaller temporary files
#----------------------------------------------------------------------------------------------------
if not skip_tagging:
    temp_files, num_records = split_pairfasta_to_temp(read1=read1, read2=read2, 
                                 temp_dir=temp_dir, 
                                 memory_limit=usr_memory_limit, threads=num_threads)
else:
    temp_files, num_records = split_fasta2_to_temp(read2=read2, 
                                 temp_dir=temp_dir, 
                                 memory_limit=usr_memory_limit, threads=num_threads)


#%% Run the barcode correction in parallel
#----------------------------------------------------------------------------------------------------
barcode_df = []  # List to collect barcode DataFrames

if not skip_tagging:

    total_files = len(temp_files)
    all_tmp_barcode_dfs = []
    temp_fastq_files = []

    # Correct Barcodes in Parallel
    with ThreadPoolExecutor(max_workers=num_threads) as executor:

        correctbc_future = {}
        for i, (temp_read1, temp_read2) in enumerate(temp_files, start=1):
            logging.info(f"Correcting barcodes for temporary file {i} out of {total_files}")
            correctbc_future[executor.submit(
                barcode_correction_pipeline,
                read2=temp_read2,
                barcode_files_dir=barcode_files_dir,
                reference=reference,
                qthreshold=qthreshold,
                num_threads=num_threads,
                usr_memory_limit=usr_memory_limit,
            )] = (temp_read1, temp_read2)

        for i, future in enumerate(as_completed(correctbc_future), start=1):
            temp_read1, temp_read2 = correctbc_future[future]
            try:
                tmp_barcode_df = future.result()
                all_tmp_barcode_dfs.append(tmp_barcode_df)  
                logging.info(f"Completed barcode correction for file {i} out of {total_files}.")
                
                # Keep track of FASTQ files for cleanup
                temp_fastq_files.append(temp_read1)
                temp_fastq_files.append(temp_read2)
            
            except Exception as e:
                logging.error(
                    f"Error correcting barcodes for read2 file stoed as: {temp_read2}\n"
                    "HELP - Use -v <verbose> for more detailed debugging"
                )
                debug_logger.debug(f"Error processing barcodes for read2 stoed as: {temp_read2}\nException: {e}")
                sys.exit(1)

    barcode_df = pd.concat(all_tmp_barcode_dfs, ignore_index=True)
    logging.info("All barcode have been corrected")

    try:
        barcode_stats, umi_stats = calculate_stats(barcode_df, expectedbc_df_list)
        logging.info(f"Stats calculation complete")
        create_heatmaps(barcode_stats, output_dir, output_base)
        create_umi_plots(umi_stats, output_dir, output_base)

    except Exception as e:
        logging.error(f"Error calculating stats and associated graphs\n"
                      "HELP - Use -v <verbose> for more detailed debugging"
                      )
        debug_logger.debug(f"Error calculating stats\nException: {e}")
        sys.exit(1)

    # Start Tagging in Parallel
    tagging_futures = []
    temp_bam_files = []

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for i, (tmp_barcode_df, (temp_read1, temp_read2)) in enumerate(zip(all_tmp_barcode_dfs, temp_files), start=1):
            temp_bam_file = tempfile.NamedTemporaryFile(delete=False, suffix=".bam", dir=temp_dir).name
            
            logging.info(f"Tagging barcodes for temporary file {i} out of {total_files}")
            tagging_futures.append(
                executor.submit(
                    tag_barcode,
                    fastq_file=temp_read1,
                    output_bam_file=temp_bam_file,
                    corrected_barcode_df=tmp_barcode_df,
                    header=bam_header
                )
            )
            temp_bam_files.append(temp_bam_file)
            
    # Wait for all tagging tasks to complete
    for i, tagging_future in enumerate(as_completed(tagging_futures), start=1):
        try:
            tagging_future.result()
            logging.info(f"Completed tagging for temporary file {i} out of {total_files}")
        except Exception as e:
            logging.error(f"Error during BAM tagging\n"
                    "HELP - Use -v <verbose> for more detailed debugging"
                )
            debug_logger.debug(f"Error during BAM tagging\nException: {e}")
            sys.exit(1)

    # Merge all temporary BAM files into the final output BAM file
    concatenate_bam_files(temp_bam_files, output_bam)
    logging.info(f"Successfully tagged all reads in {output_bam}")

    # Cleanup temporary FASTQ files
    for temp_fastq in temp_fastq_files:
        try:
            os.remove(temp_fastq)
        except Exception as e:
            logging.warning(f"Failed to delete temporary fastq file {temp_fastq}: {e}")
    debug_logger.debug("Temporary fastq files cleaned up.")

    # Cleanup temporary BAM files
    for temp_bam in temp_bam_files:
        try:
            os.remove(temp_bam)
        except Exception as e:
            logging.warning(f"Failed to delete temporary BAM file {temp_bam}: {e}")
    debug_logger.debug("Temporary BAM files cleaned up after merging.")               
    
    logging.info("Completed temporary files cleaned up")

else:
    total_files = len(temp_files)
    all_tmp_barcode_dfs = []
    temp_fastq_files = []

    # Correct Barcodes in Parallel
    with ThreadPoolExecutor(max_workers=num_threads) as executor:

        correctbc_future = {}
        for i, temp_read2 in enumerate(temp_files, start=1):
            logging.info(f"Correcting barcodes for temporary file {i} out of {total_files}.")
            correctbc_future[executor.submit(
                barcode_correction_pipeline,
                read2=temp_read2,
                barcode_files_dir=barcode_files_dir,
                reference=reference,
                qthreshold=qthreshold,
                num_threads=num_threads,
                usr_memory_limit=usr_memory_limit,
            )] = temp_read2

        for i, future in enumerate(as_completed(correctbc_future), start=1):
            temp_read2 = correctbc_future[future]
            try:
                tmp_barcode_df = future.result()
                all_tmp_barcode_dfs.append(tmp_barcode_df)  
                logging.info(f"Completed barcode correction for file {i} out of {total_files}.")
                
                # Keep track of FASTQ files for cleanup
                temp_fastq_files.append(temp_read2)
            
            except Exception as e:
                logging.error(
                    f"Error correcting barcodes for read2 file stoed as: {temp_read2}\n"
                    "HELP - Use -v <verbose> for more detailed debugging"
                )
                debug_logger.debug(f"Error processing barcodes for read2 stoed as: {temp_read2}\nException: {e}")
                sys.exit(1)

    barcode_df = pd.concat(all_tmp_barcode_dfs, ignore_index=True)
    logging.info("All barcode have been corrected.")
    try:
        barcode_stats, umi_stats = calculate_stats(barcode_df, expectedbc_df_list)
        logging.info(f"Stats calculation complete")
        create_heatmaps(barcode_stats, output_dir, output_base)
        create_umi_plots(umi_stats, output_dir, output_base)

    except Exception as e:
        logging.error(f"Error calculating stats and associated graphs\n"
                      "HELP - Use -v <verbose> for more detailed debugging"
                      )
        debug_logger.debug(f"Error calculating stats\nException: {e}")
        sys.exit(1)

    # Cleanup temporary FASTQ files
    for temp_fastq in temp_fastq_files:
        try:
            os.remove(temp_fastq)
        except Exception as e:
            logging.warning(f"Failed to delete temporary fastq file {temp_fastq}: {e}")
    debug_logger.debug("Temporary fastq files cleaned up.")

    logging.info(f"Skipping tagging as requested")
  
#%% Save stats if --stats is inputed
#----------------------------------------------------------------------------------------------------

if args.stats:

    bc_dir=os.path.join(output_dir, f"{output_base}_barcode_stats.csv")
    barcode_stats.to_csv(bc_dir, index=False)
    logging.info(f"Barcode statistics have been saved in {bc_dir}.")

    umi_dir=os.path.join(output_dir, f"{output_base}_umi_stats.csv")
    umi_stats.to_csv(umi_dir, index=False)
    logging.info(f"UMI statistics have been saved in {umi_dir}")


#%% Finileze script
#----------------------------------------------------------------------------------------------------

# Print final statistics in the log file
tagged_reads = len(barcode_df)
efficiency = (tagged_reads / num_records) * 100 if num_records > 0 else 0.0
num_cells = len(umi_stats['Cell'].unique())

logging.info(
        f"Final info:\n"
        f"\t\t\tNumber of original reads: {num_records}\n"
        f"\t\t\tNumber of taggable reads: {tagged_reads}\n"
        f"\t\t\tEfficiency: {efficiency:.2f}%\n"
        f"\t\t\tNumber of unique barcode combinations: {num_cells}\n"
    )

logging.info("BarQC completed")
logging.info("Have fun analyzing! (>ᴗ•)❀")

# Ensure all logging is outputted
for handler in logging.root.handlers[:]:
    handler.flush()
    handler.close()