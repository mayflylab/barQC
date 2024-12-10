# barQC

**barQC** is a powerful tool for processing Split-seq reads, focusing on generating detailed statistical reports and visualizations. It provides comprehensive insights into barcode distributions, including heatmaps and UMI statistics, to facilitate quality control and data analysis. Additionally, **barQC** can extract, correct, and tag barcodes for further downstream analysis.

## Table of Contents

- [Description](#description)
- [Requirements and Installation](#requirements-and-installation)
- [Usage](#usage)
  - [Arguments](#arguments)
  - [Output Files](#output-files)
  - [Additional Files](#additional-files)
- [License](#license)
- [Contact](#contact)

## Description

**barQC** enables efficient processing of Split-seq reads by:

- Extracting, correcting, and tagging barcodes from FASTQ files.
- Generating heatmaps of barcode distributions.
- Providing detailed statistics on barcodes and UMI (Unique Molecular Identifier) counts.

The tool is designed for high-performance environments, supporting multithreading and memory management to handle large datasets.

## Requirements and Installation

Clone the repository and navigate to the project directory:

```sh
git clone https://github.com/mayflylab/barQC
cd barQC
```

### Set Up the Conda Environment

Install the required dependencies using the provided Conda environment file:

```sh
conda env create -f barQC_conda_recipe.yml
conda activate barQC
```

### Install bbmap

`bbmap` is a required dependency for alignment. Unfortunately, it is not available via Conda. Follow the steps below to install it:

1. Download the latest version of BBTools from the [official BBTools website](https://sourceforge.net/projects/bbmap/).
2. Extract the downloaded archive to a directory of your choice.
3. Add the `bbmap` directory to your system’s `PATH`. For example:

    ```sh
    export PATH=/path/to/bbmap:$PATH
    ```

4. Test the installation by running:

    ```sh
    bbmap.sh
    ```

If the command runs successfully, `bbmap` is installed and ready to use.

For detailed installation instructions, refer to the [BBTools Installation Guide](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/).

## Usage

Run the tool with the following command:

```sh
python barQC.py -f1 <read1 fastq> -f2 <read2 fastq> -o <output name> [optional arguments]
```

Example:

```sh
python barQC.py -f1 data/clean_reads_1.fastq -f2 data/clean_reads_2.fastq -o results/output -b barcode_list -t 24 --stats
```

### Arguments

Required Arguments

- `-f1, --read1_fastq`: Path to the FASTQ file for read 1 (mRNA reads without adapters).
- `-f2, --read2_fastq`: Path to the FASTQ file for read 2 (barcode reads without adapters).
- `-o, --output_name`: Path and prefix for the output files.

Optional Arguments

- `-b, --bc_dir`: Directory with expected barcode files (default: current directory).
- `-q, --qval`: Quality threshold for barcode evaluation (default: 10).
- `-t, --threads`: Number of threads for parallel processing (default: 20).
- `-s, --stats`: Save statistics to output files.
- `-v, --verbose`: Enable verbose logging for debugging.
- `--skip_tagging`: Skip barcode tagging and only generate statistics.
- `--memory_limit`: Limit memory usage in GB (default: system estimate).

### Output Files

- `<output_name>.bam`: BAM file with tagged barcodes.
- `<output_name>.log`: Log file with processing details.
- `<output_name>_debug.log`: Debug log file (if `--verbose` mode is enabled).
- `<output_name>_barcode_stats.csv`: Detailed barcode statistics (if `--stats` mode is enabled).
- `<output_name>_umi_stats.csv`: UMI statistics with counts per unique barcode combination (if `--stats` mode is enabled).
- `<output_name>_heatmap_<barcode>.png`: Heatmap visualizations of barcode distributions per Split-seq plate.
- `<output_name>_UMI_counts_distribution.png`: Histogram of UMI counts.
- `<output_name>_UMI_percell_distribution.png`: Distribution of mean UMI counts per unique barcode combination.

### Additional Files

- **Barcode Files**: Expected barcode files and invariant sequences must be provided in the `barcode_list` directory. These files are crucial for the correct functioning of **barQC**.

  **Important:**  
  - The CSV files for expected barcodes must maintain a specific structure:
    - Columns: `WellPosition`, `Name`, and `Barcode`.
    - All barcodes in each file must have the same length.
  - Ensure the filenames match the expected format: `expected_barcodes_1.csv`, `expected_barcodes_2.csv`, `expected_barcodes_3.csv`.
  - Include the invariant linker sequences in `invariable-linker-sequences.fasta`.

  Any deviation from the required structure or file naming may lead to errors during processing.

- **Conda Environment File**: `barQC_conda_recipe.yml` contains the dependencies required to run barQC.
- **Test Data Directory**: Example FASTQ files are provided in the `test_data` directory for validating the installation.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

Developed by [Maria Rossello](https://github.com/m-rossello).  
For support or feedback, email Maria at [mariarossello@ub.edu](mailto:mariarossello@ub.edu).

Found a bug or have a feature request? Feel free to open an issue on the [GitHub repository](https://github.com/mayflylab/barQC).