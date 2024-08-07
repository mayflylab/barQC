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
python barQC.py -c <clean_reads.bam> -o <output.bam> [optional arguments]
```

Example:

```sh
python barQC.py -c data/clean_reads.bam -o results/output.bam -b barcode_list -q 20 -t 10 -v
```

### Arguments

- `-c, --clean_reads`: Path to the clean reads BAM file (required)
- `-o, --output_bam`: Path to the output BAM file (required)
- `-b, --bc_dir`: Directory with expected barcode files and the invariant sequence (default: barcode_list)
- `-q, --qval`: Quality threshold (default: 10)
- `-t, --threads`: Number of threads for parallel processing (default: 20)
- `-v, --verbose`: Enable verbose logging

### Output Files

- `<output.bam>`: BAM file with tagged barcodes.
- `<output_base>.log`: Log file with processing information.
- `<output_base>_debug.log`: Debug log file (if verbose mode is enabled).
- `<output_base>_stats.log`: Statistics log file with processing metrics.
- `<output_base>_heatmap_<barcode>.png`: Heatmap visualizations of barcode distribution per Split-seq plate.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

This project was created by [Maria Rossello](https://github.com/m-rossello).
For more info, support, or to send your favorite cat memes, email Maria at [mariarossello@ub.edu](mailto:mariarossello@ub.edu)
