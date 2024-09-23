# barQC

A tool for processing Split-seq reads to extract, correct, and tag barcodes. It also generates statistical reports and visualizations of the barcode distribution per Split-seq plate.

## Table of Contents

- [Description](#description)
- [Requirements and Installation](#requirements-and-installation)
- [Usage](#usage)
- [License](#license)
- [Contact](#contact)

## Description

This tool processes Split-seq reads to extract, correct, and tag barcodes. Additionally, it generates statistical reports and visualizations of the barcode distribution per Split-seq plate.

## Requirements and Installation

Clone the repository:

```sh
git clone https://github.com/mayflylab/barQC
cd barQC
```

To make sure you have all the dependencies create and activate the conda environment:

```sh
conda env create -f barQC_conda_recipe.yml
conda activate barQC
```

bbmap is used by the python script but it's intallation is not available troght conda if you don't have it already on your system please check the BBTools [Installation Guide](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/).

## Usage

To run the tool, use the following command:

```sh
python barQC.py -f1 <read1 fastq> -f2 <read2 fastq> -o <output name> [optional arguments]
```

Example:

```sh
    python barQC.py -f1 data/clean_reads_1.fastq - data/clean_reads_2.fastq -o results/output -b barcode_list -t 24 --stats
```

### Arguments

- `-f1, --read1_fastq`: Path to the fastq file of read1 (required)
- `-f2, --read2_fastq`: Path to the fastq file of read2 containing the barcode sequence (required)
- `-o, --output_name`: Path to the output files (required)
- `-b, --bc_dir`: Directory with expected barcode files (default: current directory)
- `-q, --qval`: Quality threshold (default: 10)
- `-t, --threads`: Number of threads for parallel processing (default: 20)
- `-s, --stats`: Save stats in an output file.
- `-v, --verbose`: Enable verbose logging. Usefull for debugging.
- `--skip_tagging`: Skip the tagging step and only generate statistics.

### Output Files

- `<output_name>.bam`: BAM file with tagged barcodes.
- `<output_name>.log`: Log file with processing information.
- `<output_name>_heatmap_<barcode>.png`: Heatmap visualizations of barcode distribution per Split-seq plate.
- `<output_name>_stats.log`: Statistics log file with processing metrics (if stats mode is enabled).
- `<output_name>_debug.log`: Debug log file (if verbose mode is enabled).

### Additional Files

- **Barcode Files**. The repository contains a directory (`barcode_list`) with expected barcode files and the invariant sequence used for tagging. Proper configuration of these files is crucial for the correct functioning of the barQC tool.
- **Conda Environment File**. The `barQC_conda_recipe.yml` file provides the necessary dependencies for the tool. Use this file to create and activate the conda environment required to run barQC
- **Test Data Directory**. The `test_data` directory contains example input files for testing the tool. This can be used to verify that the installation and basic functionality of the tool are correct.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

This project was created by [Maria Rossello](https://github.com/m-rossello).
For more info, support, or any feedback, email Maria at [mariarossello@ub.edu](mailto:mariarossello@ub.edu)
