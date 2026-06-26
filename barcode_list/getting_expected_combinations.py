import pandas as pd
import itertools

def generate_barcode_combinations(file1, file2, file3, output_file):
    """
    Generate all combinations of barcodes from three files and 
    save them as "barcode1barcode2barcode3".

    Each input file must have a column 'Barcode' (third column).

    Parameters
    ----------
    file1, file2, file3 : str
        Paths to CSV files containing barcodes.
    output_file : str
        Path to the output text file.
    """
    # Read the barcode column (third column)
    barcodes1 = pd.read_csv(file1)["Barcode"].tolist()
    barcodes2 = pd.read_csv(file2)["Barcode"].tolist()
    barcodes3 = pd.read_csv(file3)["Barcode"].tolist()

    # Write combinations
    with open(output_file, "w") as f:
        for b1, b2, b3 in itertools.product(barcodes1, barcodes2, barcodes3):
            f.write(f"{b1}{b2}{b3}\n")

# Example usage:
generate_barcode_combinations("expected_barcodes_1.csv", "expected_barcodes_2.csv", "expected_barcodes_3.csv", "all_combinations.txt")
