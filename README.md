# TraDb

Database and tool for the identification and annotation of plasmid transfer genes in bacterial genome sequences. The tool can process whole genome sequences, plasmids and ORFs as input. When whole genomes are submitted, the plasmid contigs are identified first, then the predicted ORFs are matched against our curated database. The database is located in the directroy `tra_db`. The file `screen_tradb.py` tries to download it automatically, but if this fails for some reason, please download all files, unzip them and save them in `$HOME/tra_db/` (default path). The tool allows you to search for and analyze sequences to find and predict genes involved in horizontal gene transfer and to track conjugable plasmids. Dependencies are the following: `prodigal 2.6.3, bedtools 2.31.0, platon 1.7, wget, mmseqs2 13.45111`.

The only mandatory and positional argument to start the script is the input, and the analysis can be started as follows: `screen_tradb.py INPUT.fasta`. The complete list of options can be found in the help menu (`screen_tradb.py -h`):
```
                                          _____          ____  _
           ___  ___ _ __ ___  ___ _ __   |_   _| __ __ _|  _ \| |__
          / __|/ __| '__/ _ \/ _ \ '_ \    | || '__/ _` | | | | '_ \
          \__ \ (__| | |  __/  __/ | | |   | || | | (_| | |_| | |_) |
          |___/\___|_|  \___|\___|_| |_|   |_||_|  \__,_|____/|_.__/ v0.1
    
usage: screen_tradb.py [-h] [--type {dna,aa,both}] [--reading_frames] [--plasmid_prediction] [--number_of_processors NUMBER_OF_PROCESSORS] [--output_directory OUTPUT_DIRECTORY] [--platon_db_dir PLATON_DB_DIR]
                       [--transfer_gene_db_dir TRANSFER_GENE_DB_DIR] [--min_aa_identity MIN_AA_IDENTITY] [--min_aa_coverage MIN_AA_COVERAGE] [--min_aa_evalue MIN_AA_EVALUE] [--min_dna_identity MIN_DNA_IDENTITY]
                       [--min_dna_coverage MIN_DNA_COVERAGE] [--min_dna_evalue MIN_DNA_EVALUE] [--version]
                       input_fasta

Screen conjugational transfer genes in plasmid sequences.

positional arguments:
  input_fasta           Input fasta file. Mandatory input.

options:
  -h, --help            show this help message and exit
  --type {dna,aa,both}, -t {dna,aa,both}
                        Type of search to be run.
  --reading_frames, -ORF
                        If the input consits of ORFs. Use this only if you alredy predicted the ORFs of the input. If set, you should also consider to change --type appropriately. Default is False.
  --plasmid_prediction, -p
                        If plasmid prediction of the input should be attempted. The output of platon containing the plasmids will be screened for transfer genes
  --number_of_processors NUMBER_OF_PROCESSORS, -np NUMBER_OF_PROCESSORS
                        Number of processes invoked for the search.
  --output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
                        Output directory to store results. Output files will be named using the basename of the input file.
  --platon_db_dir PLATON_DB_DIR, -pdb PLATON_DB_DIR
                        Where the reference database of platon can be found. Default is $HOME/platon_db. If it can not be found automatic install will be attempted.
  --transfer_gene_db_dir TRANSFER_GENE_DB_DIR, -tdb TRANSFER_GENE_DB_DIR
                        Where the reference conjugational transfer gene database can be found. Default is $HOME/tra_db. If it can not be found automatic install will be attempted.
  --min_aa_identity MIN_AA_IDENTITY, -mai MIN_AA_IDENTITY
                        Minimum percent of identity when screening amino acids (0-100). Default is 80.
  --min_aa_coverage MIN_AA_COVERAGE, -mac MIN_AA_COVERAGE
                        Minimum percent of identity when screening amino acids (0-100). Default is 80.
  --min_aa_evalue MIN_AA_EVALUE, -mae MIN_AA_EVALUE
                        Minimum e-value when screening amino acids. Default is 1e-50.
  --min_dna_identity MIN_DNA_IDENTITY, -mdi MIN_DNA_IDENTITY
                        Minimum percent of identity when screening DNA sequences (0-100). Default is 80.
  --min_dna_coverage MIN_DNA_COVERAGE, -mdc MIN_DNA_COVERAGE
                        Minimum percent of identity when screening DNA sequences (0-100). Default is 80.
  --min_dna_evalue MIN_DNA_EVALUE, -mde MIN_DNA_EVALUE
                        Minimum e-value when screening DNA sequences. Default is 1e-50.
  --version, -V         Output version number and exit
```

THIS IS A WORK IN PROGRESS. USE AT YOUR OWN RISK AND DOUBLE CHECK YOUR RESULTS.
