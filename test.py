#!/usr/bin/env python

import subprocess
import filecmp
import os
import difflib

# Command to run
# Loop through all .fasta files in the ./test/ directory
fasta_files = [f for f in os.listdir('./test/') if f.endswith('.fasta')]
for fasta_file in fasta_files:
    command = f'./screen_tradb.py -np 4 ./test/{fasta_file} -t both'
    print(f"Running command: {command}")
    # Run the command
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command for {fasta_file}: {e}")
        exit(1)
print(f"Current working directory: {os.getcwd()}")
# List the contents of the current working directory
print("Contents of the working directory:")
for item in os.listdir(os.getcwd()):
    print(f" - {item}")
# Run the command
try:
    subprocess.run(command, check=True, shell=True)
except subprocess.CalledProcessError as e:
    print(f"Error running command: {e}")
    exit(1)

# Directory containing expected output
expected_output_dir = 'test/test_out'

# Directory containing actual output (assuming it is generated in the current directory)
actual_output_dir = './tra_out'

# Function to sort a TSV file by the first column
def sort_tsv(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    header = lines[0] if lines and '\t' in lines[0] else None
    data_lines = lines[1:] if header else lines
    sorted_lines = sorted(data_lines, key=lambda x: x.split('\t')[0])
    return [header] + sorted_lines if header else sorted_lines

# Compare files
success = True
for filename in os.listdir(expected_output_dir):
    # Skip .log files
    if filename.endswith('.log') or filename.endswith('raw_hits.tsv'):
        continue

    expected_file = os.path.join(expected_output_dir, filename)
    actual_file = os.path.join(actual_output_dir, filename)

    if not os.path.exists(actual_file):
        print(f"Missing output file: {filename}")
        success = False
        continue

    # Sort and compare TSV files
    if filename.endswith('.tsv'):
        expected_sorted = sort_tsv(expected_file)
        actual_sorted = sort_tsv(actual_file)
        if expected_sorted != actual_sorted:
            success = False
            print(f"Differences found in file: {filename}")
            diff = difflib.unified_diff(expected_sorted, actual_sorted, 
                                        fromfile='expected', tofile='actual')
            print(''.join(diff))
    else:
        if not filecmp.cmp(expected_file, actual_file, shallow=False):
            success = False
            print(f"Differences found in file: {filename}")
            with open(expected_file, 'r') as ef, open(actual_file, 'r') as af:
                expected_lines = ef.readlines()
                actual_lines = af.readlines()
                diff = difflib.unified_diff(expected_lines, actual_lines, 
                                            fromfile='expected', tofile='actual')
                print(''.join(diff))

if success:
    print("Test successful: All outputs match expected results.")
else:
    print("Test failed: Differences found.")