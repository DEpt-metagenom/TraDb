#!/usr/bin/env python

import os
import subprocess
import sys

home_dir = os.path.expanduser('~')

output_dir = sys.argv[1] if len(sys.argv) > 1 else 'home_dir'

os.makedirs(output_dir, exist_ok=True)

platon_db_dir = os.path.join(output_dir, 'platon_db')
tra_db_dir = os.path.join(output_dir, 'tra_db')


def run_shell_cmd(cmd):
    """
    Run a shell command and handle errors.

    Args:
        cmd (str): The shell command to run.

    Returns:
        None
    """
    try:
        print(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running command: {cmd}")
        print(f"Error message: {e}")
        raise
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise


def download_and_extract(url, output_dir, output_file, arch_type='tar'):
    """
    Download and extract a file from a URL.

    Args:
        url (str): The URL to download the file from.
        output_dir (str): The directory to extract the file into.
        output_file (str): The name of the output file to save the downloaded content.
        arch_type (str, optional): The type of archive to extract. Defaults to 'tar'.
                                   Supported types are 'tar', 'gz', and None.

    Returns:
        None

    Outputs:
        - The downloaded file is saved to the specified output directory.
        - The file is extracted into the output directory if it is an archive.
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    wget_cmd = f'wget {url} -O {output_file}'
    if arch_type == 'tar':
        arch_cmd = f'mkdir -p {output_dir} && tar -xzf {output_file} -C {output_dir}'
    elif arch_type == 'gz':
        arch_cmd = f'mkdir -p {output_dir} && gunzip -f {output_file}'
    elif arch_type is None:
        arch_cmd = f'echo "metadata database downloaded"'
    run_shell_cmd(wget_cmd)
    run_shell_cmd(arch_cmd)


def download_all_databases(tra_db_dir, platon_db_dir):
    """
    Download all required databases for the script.

    Args:
        tra_db_dir (str): Directory where the TraDb database will be stored.
        platon_db_dir (str): Directory where the Platon database will be stored.

    Returns:
        None
    """
    print("Downloading all required databases...")

    # Download TraDb DNA database
    download_and_extract(
        'https://zenodo.org/records/15539975/files/Tra_db_sequences_dna.fa.gz',
        tra_db_dir,
        f'{tra_db_dir}/Tra_db_sequences_dna.fa.gz',
        arch_type='gz'
    )

    # Download TraDb AA database
    download_and_extract(
        'https://zenodo.org/records/15539975/files/Tra_db_sequences_aa.fa.gz',
        tra_db_dir,
        f'{tra_db_dir}/Tra_db_sequences_aa.fa.gz',
        arch_type='gz'
    )

    # Download TraDb metadata for DNA
    download_and_extract(
        'https://zenodo.org/records/15539975/files/Tra_db_taxa_positons_dna.db',
        tra_db_dir,
        f'{tra_db_dir}/Tra_db_taxa_positons_dna.db',
        arch_type=None
    )

    # Download TraDb metadata for AA
    download_and_extract(
        'https://zenodo.org/records/15539975/files/Tra_db_taxa_positons_aa.db',
        tra_db_dir,
        f'{tra_db_dir}/Tra_db_taxa_positons_aa.db',
        arch_type=None
    )

    # Download Platon database
    download_and_extract(
        'https://zenodo.org/record/4066768/files/db.tar.gz',
        platon_db_dir,
        f'{platon_db_dir}/db.tar.gz',
        arch_type='tar'
    )

    print("All databases have been downloaded successfully.")

def main():
    """
    Main function to execute the script.

    Returns:
        None
    """
    print("Starting database download process...")
    download_all_databases(tra_db_dir, platon_db_dir)
    print("Database download process completed successfully.")

if __name__ == '__main__':
    try:
        main()
    except Exception as e:  # Catch any exceptions that occur during the execution
        print(f"An error occurred: {e}")
        sys.exit(1)