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

# Parse the arguments
args = parser.parse_args()


#%%% Check all the requirements are provided
#--------------------------------------------------------------------------------------------------------

# Check if all required arguments are provided and validate paths
required_args_path = {'read1_fastq': args.read1_fastq, 'read2_fastq': args.read2_fastq,'bc_dir': args.bc_dir}
missing_args_path = [arg for arg, path in required_args_path.items() if (path is None or not os.path.isfile(path)) and arg != 'bc_dir']
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
required_args = ['output_name']
missing_args = [arg for arg in required_args if getattr(args, arg) is None]
if missing_args:
    parser.print_help()
    print(f"\nError: Missing required argument: {', '.join(missing_args)}\n")
    exit(1)


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


#%% FUNCTIONS
#########################################################################################################


#%%% Auxiliary functions to map the barcodes
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
        if debug_logger:
            debug_logger.debug(f"Impossible to correct the barcode in {qname} for {bc_type}")
    
    return None


#%%% Functions to open files
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

        return string_fastq, df_fastq

    except Exception as e:
        logging.error(
            f"Failed to parse fastq file {fastq_file}: {e}. "
            "Ensure that the read2 fastq file is correctly formatted and accessible. "
            "You may also want to check if the file is corrupted or improperly indexed."
        )
        sys.exit(1)

#%%% Functions to obtain the align barcodes
#--------------------------------------------------------------------------------------------------------

def run_bbmap(string_fastq, reference, cores):

    """
    Runs bbmap.sh to align sequences and find barcodes.

    Args:
        fastq_as_df (pd.DataFrame): DataFrame containing BAM data.
        reference (str): Path to the reference genome.
        cores (int): Number of threads to use for parallel processing.

    Returns:
        bytes: Mapped BAM data in SAM format.
    """
    logging.info("Aligning sequence to find barcodes. This can take a while.")

    try:
        # Join the FASTQ entries into a single string and encode it
        fastq_data = '\n'.join(string_fastq).encode()

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
            text=False  # Important for binary input/output
        )

        bbmap_stdout, bbmap_stderr = process.communicate(input=fastq_data)

        logging.info("Sequence aligned to find barcodes")
            
        process.wait()

        return bbmap_stdout
    
    except subprocess.CalledProcessError as e:
        logging.error(
            f"Error occurred while running bbmap.sh: {e}. "
            "This may be due to issues with the reference genome, input FASTQ data, or system resources. "
            "Check the stderr output for more details:\n{e.stderr.decode()}"
        )
        logging.error(
            "Ensure that the reference file exists, is correctly formatted, "
            "and that there is sufficient memory and CPU available."
        )
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

    # print(lines)

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
                        

#%%% Process barcodes
#--------------------------------------------------------------------------------------------------------

from tqdm import tqdm

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

    # Process in chunks
    num_chunks = int(np.ceil(len(clean_df) / chunk_size))

    with tqdm(total=len(clean_df), desc="Processing barcodes") as pbar:
        for i in range(num_chunks):
            chunk = clean_df[i * chunk_size:(i + 1) * chunk_size]

            for entry in chunk.itertuples(index=False):
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

            pbar.update(len(chunk))  # Update progress bar after processing each chunk

    # Convert the list of barcodes and statistics to DataFrames
    barcode_df = pd.DataFrame(barcode_ls)

    if barcode_stats_list:
        barcode_stats = pd.DataFrame(barcode_stats_list)
        barcode_stats = barcode_stats.groupby(['WellPosition', 'Name', 'Barcode'], as_index=False).sum()
    else:
        barcode_stats = pd.DataFrame(columns=['WellPosition', 'Name', 'Barcode', 'Count'])

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
    
    # Save the heatmap to a file
    plt.savefig(filename)
    plt.close()


#%%% Tag final bam
#--------------------------------------------------------------------------------------------------------

def tag_barcode(fastq_file, output_bam_file, corrected_barcode_df, chunk_size=10000):
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
                        a = pysam.AlignedSegment()
                        a.query_name = entry.name
                        a.query_sequence = entry.sequence
                        a.query_qualities = pysam.qualitystring_to_array(entry.quality)
                        a.flag = 4  # Unmapped flag

                        # Tag barcodes if available
                        if entry.name in barcode_dict:
                            barcode_info = barcode_dict[entry.name]
                            xd, xe, xf, xm = barcode_info['BC1'], barcode_info['BC2'], barcode_info['BC3'], barcode_info['UMI']
                            xc = xd + xe + xf

                            # Add tags to the read
                            a.set_tag('XD', xd, value_type='Z')
                            a.set_tag('XE', xe, value_type='Z')
                            a.set_tag('XF', xf, value_type='Z')
                            a.set_tag('XC', xc, value_type='Z')
                            a.set_tag('XM', xm, value_type='Z')

                        # Add the read to the chunk
                        chunk.append(a)

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

        logging.info(f"Successfully tagged to {output_bam_file}")

    except Exception as e:
        logging.error(
            f"Failed to tag barcodes: {e}\n"
            "HELP - This might be due to an issue with read 1 fasta or memory limitations. "
            "Please check read1 fasta and system memory limitations.\n"
            "HELP - Use -v <verbose> for more detailed debugging"
        )


#%% PROCESS DATA
#########################################################################################################

def main():

    # Define the variables
    read1 = args.read1_fastq
    read2 = args.read2_fastq
    num_threads = args.threads
    barcode_files_dir = args.bc_dir
    qthreshold = args.qval
    output_bam = os.path.join(output_dir, f"{output_base}.bam")
    reference = os.path.join(barcode_files_dir, "invariable-linker-sequences.fasta")
    skip_tagging = args.skip_tagging  

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

    if clean_df is not None:
        
        # Use ThreadPoolExecutor for running bbmap and obtaining coordinates
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            future_bbmap = executor.submit(run_bbmap, string_fastq, reference, num_threads)

            try:
                mapped_bam_data = future_bbmap.result()
                logging.info("Barcode alignment to invariant sequence completed successfully.")

                # Calculate barcode coordinates
                future_map_dict = executor.submit(parse_and_compute_coordinates, mapped_bam_data)
                map_dict = future_map_dict.result()
                logging.info("Coordinates computed")
                gc.collect()

            except Exception as e:
                logging.error(
                    f"Failed to align barcodes or parse coordinates: {e}.\n"
                    "HELP - Possible causes include issues with the reference file, improper barcode extraction, or memory limitations. "
                    "Ensure all input files are correct and check system memory limitations.\n"
                    "HELP - Use -v <verbose> for more detailed debugging"
                )
                sys.exit(1)

    original_reads = len(clean_df) if clean_df is not None else 0
    tagged_reads = 0
    efficiency = 0.0
    num_cells = 0

    if map_dict and not clean_df.empty:

        # Use ThreadPoolExecutor for barcode processing and calculating stats
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

            except Exception as e:
                logging.error(
                    f"Failed to process barcodes: {e}\n"
                    "HELP - This might be due to an issue with the input data, such as not finding the barcodes. "
                    "Please check the input files and the barcode directory\n"
                    "HELP - Use -v <verbose> for more detailed debugging"
                )
                sys.exit(1)

    if not skip_tagging:

        # Use ThreadPoolExecutor to tag read1 with barcodes and write a bam
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            executor.submit(tag_barcode, fastq_file=read1,output_bam_file=output_bam,corrected_barcode_df=barcode_df)

    else:
        logging.info("Skipping the tagging step as requested.")

    

    # Print stats_df to stats_logger
    if stats_logger:
        stats_df = stats_df.sort_values(by=['Name', 'WellPosition'])
        stats_df_dropped = stats_df.drop(columns=['Column', 'Row'])
        stats_logger.info(stats_df_dropped.to_string(index=False))

    logging.info("Barcode processing pipeline completed.")
    logging.info("Have fun analyzing! (>ᴗ•)❀\n")

    # Print final statistics
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