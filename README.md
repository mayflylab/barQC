# barQC

A tool for processing Split-seq reads to extract, correct, and tag barcodes. It also generates statistical reports and visualizations of the barcode distribution per Split-seq plate.

## Table of Contents

- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Arguments](#arguments)
- [Output Files](#output-files)
- [Contributing](#contributing)
- [License](#license)
- [Author](#author)

## Description

This tool processes Split-seq reads to extract, correct, and tag barcodes. The reads must be clean of Illumina adapters. Additionally, it generates statistical reports and visualizations of the barcode distribution per Split-seq plate.

## Installation

Clone the repository:

```sh
git clone https://github.com/yourusername/barQC.git
cd barQC
```

## Usage

To run the tool, use the following command:

```sh
python barQC.py -c <clean_reads.bam> -o <output.bam> [optional arguments]
```

### Example:

```sh
python barQC.py -c data/clean_reads.bam -o results/output.bam -b data/barcode_list -q 20 -t 10 -v
```

## Arguments

- `-c, --clean_reads`: Path to the clean reads BAM file (required)
- `-o, --output_bam`: Path to the output BAM file (required)
- `-b, --bc_dir`: Directory with expected barcode files (default: current directory)
- `-q, --qval`: Quality threshold (default: 10)
- `-t, --threads`: Number of threads for parallel processing (default: 20)
- `-v, --verbose`: Enable verbose logging

## Output Files

- `<output.bam>`: BAM file with tagged barcodes.
- `<output_base>.log`: Log file with processing information.
- `<output_base>_debug.log`: Debug log file (if verbose mode is enabled).
- `<output_base>_stats.log`: Statistics log file with processing metrics.
- `<output_base>_heatmap_<barcode>.png`: Heatmap visualizations of barcode distribution per Split-seq plate.

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository.
2. Create a new branch: `git checkout -b my-feature-branch`
3. Make your changes and commit them: `git commit -m 'Add new feature'`
4. Push to the branch: `git push origin my-feature-branch`
5. Submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

This project was created by [Maria Rossello](https://github.com/m-rossello).
