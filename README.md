# TraDb

TraDb is a database and toolset for the identification and annotation of plasmid transfer genes in bacterial genome sequences. It supports whole genome sequences, plasmids and ORFs as input. TraDb allows searching and analysing sequences to predict genes involved in horizontal gene transfer and to track conjugative plasmids.

## Installation

The main script for searching for plasmid transfer genes is `screen_tradb.py`. The hard dependencies are the following:
- `prodigal 2.6.3`
- `bedtools 2.31.0`
- `platon 1.7`
- `wget`
- `mmseqs2 15.6f452`. 

To install TraDb and its dependencies, clone this repository using git and navigate to the TraDb directory:
```
git clone https://github.com/DEpt-metagenom/TraDb.git
cd TraDb
```

### Conda

You can create the environment using either of the following:
```
conda env create --file conda_env.yml
```

Or

```
conda create -n tradb --file conda_env_explicit.txt
```

If you prefer to install dependencies in the currently active environment you should run:
```
conda install --file conda_env_explicit.txt
```

To make `screen_tradb.py` globally available, it should be added to `$PATH`, e.g. assuming you are in the directory where `screen_tradb.py` is located, you should run:
```
export PATH=$(realpath ./):$PATH
```

Optionally, the export can be made permanent by adding the above line to your `.bashrc`, which is probably located in `$HOME/.bashrc`.

### Databases

The tool relies on multiple databases:
1. [transfer gene database] (https://zenodo.org/records/15539975), which contains DNA and AA reference sequences with metadata and was developed specifically for this tool.

2. The [platon database](https://zenodo.org/records/4066768) for the identification of plasmids in whole genome sequences.

On the first run, screen_tradb.py attempts to download all required databases for the run. Alternatively, you can also run:

```
screen_tradb.py --download-database
```

By default, the databases are saved in $HOME/tra_db. This can be changed using `--transfer_gene_db_dir` and `--platon_db_dir`, which affects both the download of the data and the location where `screen_tradb.py` looks for the reference.

You can also download the required databases with `download_all_db.py`. By default, it stores the downloaded files in `$HOME/tra_db/{trad_db,platon}`. The location of the downloaded files can be specified as a positional argument, e.g. `download_all_db.py path/to/my_local_copy_of_tra_db`.

### Test

The installation can be tested by `./test.py`˛The test script takes three plasmid sequences from the `test` directory (`AP000342.1.fasta`, `AP001918.1.fasta`, `BN000925.1.fasta`) and executes `screen_tradb.py` for all of them. Then it compares the results with the expected output files in `test/test_out/`. If the results match, the script should report: "`Test successful: All outputs match the expected results`". In case of differences, the test fails and the dependencies should be double checked. The test script reports differences and `screen_tradb.py` reports missing dependencies. If the test script reports "`Test failed: Differences found.`", please check the test log carefully.

### Containers

This tool is also distributed as an image hosted on DockerHub with all dependencies pre-installed. It can be run using Docker or Apptainer (formerly Singularity). The reference databases are not included in the image and should be downloaded prior to starting the container.

Assuming the reference database (`tra_db`) and the plasmid to be analysed are in the current directory, `screen_tradb.py` can be started using the Docker container in the following way:

```
# downloads databases to the current working directory, run only once
download_all_db.py ./ 

# starts screen_tradb.py using Docker
# input.fasta is a plasmid sequence in the current working directory
docker run -i -v $(pwd):/data/ -t deptmetagenom/tra_db:0.2.5 screen_tradb.py --transfer_gene_db_dir /data/tra_db/ /data/input.fasta
```

The Docker container can also be built locally by running `sudo docker build -t tra_db .`, provided that the Dockerfile is located in the current directory.

If you are using an HPC system or Docker is not available for you, Apptainer might be more practical to run the container. By default, Apptainer uses automount and should find the database in its default location (or at the specified path). `screen_tradb.py` can be started with Apptainer in the following way:

```
apptainer exec docker://deptmetagenom/tra_db:0.2.5 screen_tradb.py input.fasta
```

The above command pulls the image from DockerHub. If you want to save a local copy of the container as a portable `sif` file, you can execute this command:

```
apptainer pull docker://deptmetagenom/tra_db:0.2.5
```

Then `screen_tradb.py` can be started by running:

```
apptainer exec /path/to/tra_db_0.2.5.sif screen_tradb.py 
```

Alternatively, you can also build the `.sif` container locally. This gives you the opportunity to customise the container to your own requirements. To build the container locally, execute:

```
apptainer build tradb.sif apptainer.def 
```

Then the container should run as if it were pulled from DockerHub.

`.sif` files are portable and can be copied to different computers, so the tools should run identically on a local and a remote computer. If you do this, please make sure that the reference database is also copied or downloaded to the remote computer.

## Usage

The only mandatory and positional argument for starting the script is the input, assumed the input consists of plasmid sequences.

The script can be executed in three modes.

### Plasmids

This is the default mode for running `screen_tradb.py`, which does not require any options other than the defaults. In this mode it is assumed that the input plasmid sequences are in fasta format. The script predicts ORFs, which are matched against the reference database in DNA and AA format by default. The metadata (species, position of the genes in the reference genome) is then queried from an sqlite3 database. The default identity and coverage of the reference sequence queries is 80 %. The output is stored in `tra_out/`. This type of analysis (using one thread by default) can be started as follows:

`screen_tradb.py input.fasta`.

### Whole genome sequences

If whole genomes are used as input, `screen_tradb.py` identifies the plasmid contigs with `platon` before predicting ORFs and searching for transfer genes. The platon database is only required if you use this feature. This feature can be activated using the option `--plasmid_prediction`, e.g:

`screen_tradb.py --plasmid_prediction whole_genome_with_chromosome_and_plasmid.fasta`.

### ORFs

If you have already predicted the ORFs of the contigs of interest, you can provide the ORF sequences directly to `screen_tradb.py`. The feature can be switched on using the `--reading_frames` option. In this way, only the similarity search of ORFs and the query of metadata is performed. ORFs can be specified as either DNA or AA and `--type` should be selected accordingly. The option `--type auto` is only valid in this mode and will attempt to determine whether the input is DNA or AA. The identification of transfer genes in predicted ORFs can be performed as follows:

`screen_tradb.py --reading_frames --type {dna,aa,auto} ORFs.fasta`.

## Options
If you need to change the output directory, the number of threads, the minimum query identity and coverage, or the type of search, please refer to the list of options to find the appropriate option, which can be found in the help menu (by running `screen_tradb.py -h`):

```
                                          _____          ____  _
           ___  ___ _ __ ___  ___ _ __   |_   _| __ __ _|  _ \| |__
          / __|/ __| '__/ _ \/ _ \ '_ \    | || '__/ _` | | | | '_ \
          \__ \ (__| | |  __/  __/ | | |   | || | | (_| | |_| | |_) |
          |___/\___|_|  \___|\___|_| |_|   |_||_|  \__,_|____/|_.__/ v0.2.5
    
usage: screen_tradb.py [-h] [--download-database] [--type {dna,aa,both,auto}]
                       [--reading_frames] [--plasmid_prediction]
                       [--number_of_processors NUMBER_OF_PROCESSORS]
                       [--output_directory OUTPUT_DIRECTORY]
                       [--platon_db_dir PLATON_DB_DIR]
                       [--transfer_gene_db_dir TRANSFER_GENE_DB_DIR]
                       [--min_aa_identity MIN_AA_IDENTITY]
                       [--min_aa_coverage MIN_AA_COVERAGE]
                       [--min_aa_evalue MIN_AA_EVALUE]
                       [--min_dna_identity MIN_DNA_IDENTITY]
                       [--min_dna_coverage MIN_DNA_COVERAGE]
                       [--min_dna_evalue MIN_DNA_EVALUE] [--version]
                       [input_fasta]

Screening of conjugation transfer genes in plasmid sequences.

positional arguments:
  input_fasta           Input fasta file. Mandatory input unless --download-
                        database is used.

optional arguments:
  -h, --help            show this help message and exit
  --download-database, -dd
                        Download all required databases and exit.
  --type {dna,aa,both,auto}, -t {dna,aa,both,auto}
                        Type of search to be run. Choices are 'dna', 'aa',
                        'both' and 'auto'. Default is 'both'. Auto works only
                        if the input are ORFs (-ORF or --reading-frames).
  --reading_frames, -ORF
                        If the input consists of ORFs. Use this only if you
                        already predicted the ORFs of the input. If set, you
                        should also consider changing --type appropriately.
                        Default is False.
  --plasmid_prediction, -p
                        If plasmid prediction of the input should be
                        attempted. The output of platon containing the
                        plasmids will be screened for transfer genes
  --number_of_processors NUMBER_OF_PROCESSORS, -np NUMBER_OF_PROCESSORS
                        Number of processes invoked for the search.
  --output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
                        Output directory to store results. Output files will
                        be named using the basename of the input file.
  --platon_db_dir PLATON_DB_DIR, -pdb PLATON_DB_DIR
                        Where the reference database of platon can be found.
                        Default is $HOME/platon_db. If it cannot be found,
                        automatic install will be attempted.
  --transfer_gene_db_dir TRANSFER_GENE_DB_DIR, -tdb TRANSFER_GENE_DB_DIR
                        Where the reference conjugational transfer gene
                        database can be found. Default is $HOME/tra_db. If it
                        cannot be found, automatic install will be attempted.
  --min_aa_identity MIN_AA_IDENTITY, -mai MIN_AA_IDENTITY
                        Minimum percent of identity when screening amino acids
                        (0-100). Default is 80.
  --min_aa_coverage MIN_AA_COVERAGE, -mac MIN_AA_COVERAGE
                        Minimum percent of coverage when screening amino acids
                        (0-1). Default is 0.8.
  --min_aa_evalue MIN_AA_EVALUE, -mae MIN_AA_EVALUE
                        Minimum e-value when screening amino acids. Default is
                        1e-50.
  --min_dna_identity MIN_DNA_IDENTITY, -mdi MIN_DNA_IDENTITY
                        Minimum percent of identity when screening DNA
                        sequences (0-100). Default is 80.
  --min_dna_coverage MIN_DNA_COVERAGE, -mdc MIN_DNA_COVERAGE
                        Minimum percent of coverage when screening DNA
                        sequences (0-1). Default is 0.8.
  --min_dna_evalue MIN_DNA_EVALUE, -mde MIN_DNA_EVALUE
                        Minimum e-value when screening DNA sequences. Default
                        is 1e-50.
  --version, -V         Print version number and exit
```


## Output
The output of `screen_tradb.py` contains the sequences of the transfer genes, the genomic coordinates of the transfer genes and the metadata of the reference sequences (species, position of the gene in the original plasmid). The output files (if plasmids were used as input) are the following:

- `<basename>_dna.fa`, `<basename>_aa.fa`: Predicted ORFs in DNA and AA format.
- `<basename>_dna_raw_hits.tsv`, `<basename>_aa_raw_hits.tsv`: The results of mmseqs2 in tabular format using DNA and AA as input without any filters applied.
- `<basename>_dna_filtered_hits.tsv`, `<basename>_aa_filtered_hits.tsv`: The coordinates of the transfer genes with their identity and coverage values. The ID of the hits (Cluster column) indicates the gene and the similarity cluster identified when building the reference database, e.g. TraA-Cluster10 means that TraA was found and the gene belongs to similarity cluster 10. The cluster IDs are not meaningful on their own, but are used to retrieve the metadata of the hits. The query ID displays the ID of the sequence and the genomic coordinates of the gene in the following format: `contig_ID:start-end`.
- `<basename>_dna_hit_sequences.fa`, `<basename>_aa_hit_sequences.fa`: The sequence of the identified transfer genes in DNA and AA format.
- `<basename>_dna_hit_data.tsv`, `<basename>_aa_hit_data.tsv`: The metadata of the identified transfer genes. The filtered list of hits is supplemented by the list of taxa and the position of the transfer genes in the original plasmid sequence.
- `<basename>_dna_mmseqs.log`, `<basename>_aa_mmseqs.log`: MMSeqs2 log of the similarity search using DNA and AA.
- `<basename>.bed`: Coordinates of the ORFs in the input sequence.
- `<basename>.gff`: The output of Prodigal (used to identify ORFs) in gff format.

## Citation

`screen_tradb.py` relies on multiple external tools to find transfer genes. If you use this tool, please do not forget to cite the original work:

Hyatt, D., Chen, G. L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics, 11, 1-11.

Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842.

Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology, 35(11), 1026-1028.

If you used platon to identify plasmids, please cite:

Schwengers, O., Barth, P., Falgenhauer, L., Hain, T., Chakraborty, T., & Goesmann, A. (2020). Platon: identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein sequence-based replicon distribution scores. Microbial genomics, 6(10), e000398.
