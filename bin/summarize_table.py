import os
import glob
import pandas as pd
import argparse

def concatenate_files(input_dir, output_file, prefix):
    all_files = glob.glob(os.path.join(input_dir, "*"))
    li = []

    for filename in all_files:
        base_name = os.path.basename(filename)
        modified_name = base_name.replace(args.prefix, '')
        
        df = pd.read_csv(filename, index_col=None, header=0, sep='\t')
        df['ID'] = modified_name
        li.append(df)

    frame = pd.concat(li, axis=0, ignore_index=True)
    frame = frame[['ID'] + [col for col in frame.columns if col != 'ID']]
    frame.to_csv(output_file, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Concatenate multiple TSV files into one, adding ID (from filename) as a column.')
    parser.add_argument('-i', '--input-dir', required=True, help='Directory containing the TSV files to concatenate.')
    parser.add_argument('-p', '--prefix', required=True, help='Prefix to be removed from the filenames.')
    parser.add_argument('-o', '--output-file', required=True, help='Path to the output file.')
    args = parser.parse_args()

    concatenate_files(args.input_dir, args.output_file, args.prefix)