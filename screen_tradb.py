#!/usr/bin/env python

import argparse
import pandas as pd
import os
import re
import sys

version_number = str("screen_tradb version 0.1")

home_dir = os.path.expanduser('~')

def run_platon(infile, basename, outdir, num_proc, platon_db_dir):
    """define function to run platon"""

    platon_cmd = f'platon --db {home_dir}/platon_db/db/ --output {basename}_platon --mode accuracy --threads {num_proc} {infile}'

    if not os.path.exists(platon_db_dir):
        wget_cmd = f'wget https://zenodo.org/record/4066768/files/db.tar.gz -O {home_dir}/platon_db.tar.gz'
        tar_cmd = f'mkdir {home_dir}/platon_db && tar -xzf {home_dir}/platon_db.tar.gz -C {home_dir}/platon_db'
        os.system(wget_cmd)
        os.system(tar_cmd)

    print(f'predicting plasmids with platon using the command: {platon_cmd}')
    os.system(platon_cmd)

def run_prodigal(infile, basename, outdir):
    """define function to run prodigal"""

    prodigal_cmd = f'prodigal -i {infile} -p meta -f gff -a {outdir}/{basename}_aa.fa -o {outdir}/{basename}.gff 2> /dev/null'
    os.system(prodigal_cmd)
    print("ORF prediction complete")

    # could parse as fasta but as the replaced characters must not appear in AA sequences this solution will do for now
    # TODO remove partial!
    pattern = r'>([\w\.]+)_(\d+) # (\d+) # (\d+).+'

    with open(f'{outdir}/{basename}_aa.fa', 'r') as file:
        amino_acids = file.read()
        amino_acids = re.sub(pattern, r'>\1:\3-\4', amino_acids)

    with open(f'{outdir}/{basename}_aa.fa', 'w') as file:
        file.write(amino_acids)

    # sometimes prodigal fails to output the DNA (no problem with AA), below there is a workaround to this rare problem
    filtered_lines = []
    with open(f'{outdir}/{basename}.gff', 'r') as gff_file:
        for line in gff_file:
            if not line.startswith("#") and "partial=01" not in line:
                fields = line.strip().split("\t")
                start = str(int(fields[3]) - 1)  # correct 1-based coordinates
                filtered_lines.append((fields[0], start, fields[4]))

    with open(f'{outdir}/{basename}.bed', 'w') as bed_file:
        for entry in filtered_lines:
            bed_file.write("\t".join(entry) + "\n")

    bedtools_cmd = f'bedtools getfasta -fi {infile} -bed {outdir}/{basename}.bed > {outdir}/{basename}_dna.fa'
    os.system(bedtools_cmd)

def filter_hits(raw_hits, filtered_hits):
    """define function to filter blast hits by score"""

    hits = pd.read_csv(raw_hits, sep='\t')

    hits['score'] = pd.to_numeric(hits['score'], errors='coerce')

    filtered = hits.groupby('hit', as_index=False).apply(lambda group: group.loc[group['score'].idxmax()])

    filtered.reset_index(drop=True, inplace=True)

    filtered.to_csv(filtered_hits, sep='\t', index=False)

# TODO frequency of reference taxa / gene
def get_metadata(filtered_hits, tra_amino_data, hit_data):
    """function to get metadata of hits, i.e. taxa the hit can be associated with and the position of reference genes(contig:start-end) using 0-based coordinates"""

    hits = pd.read_csv(filtered_hits, sep='\t')
    metadata_columns = ['cluster', 'taxa', 'position']
    metadata = pd.read_csv(tra_amino_data, sep='\t', names=metadata_columns)

    # for testing only
    # matching_values = metadata['cluster'].isin(hits['cluster'])
    # result = metadata[matching_values]
    # print(result)

    merged_data = pd.merge(hits, metadata, on='cluster')
    merged_data.to_csv(hit_data, sep="\t", index=False)

def filter_fasta(orf_fasta, hit_data, out_fasta):
    """function to subset ORFs in fasta format with Tra hits"""

    sequences = []
    seq = None
    name = ''
    filtered = []
    hits = pd.read_csv(hit_data, sep='\t')

    tra_ids = hits['hit'].tolist()

    with open(orf_fasta, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if seq is not None:
                    sequences.append((name, seq))
                name = line[1:]
                seq = ''
            else:
                seq += line
        if seq is not None:
            sequences.append((name, seq))

    for name, seq in sequences:
        if name in tra_ids:
            tra_row = hits[hits['hit'] == name]
            tra_id = tra_row.loc[:, 'cluster'].values
            formatted_tra_id = ', '.join(tra_id)
            filtered.append((f'{formatted_tra_id} {name}', seq))

    with open(out_fasta, 'w') as file:
        for name, sequence in filtered:
            file.write(f'>{name}\n{sequence}\n')

def run_blastp(amino_acids_db, tra_db_dir, tra_amino_acids, tra_amino_data, basename, outdir, num_proc, mai, mac, mae):
    """define function to run blastp"""

    if not os.path.exists(f'{amino_acids_db}.ndb'):
        make_amino_acids_db_cmd = f'makeblastdb -in {outdir}/{basename}_aa.fa -dbtype prot -out {amino_acids_db} 1> /dev/null'
        os.system(make_amino_acids_db_cmd)
        print("Database of ORFs in amino acid format created")

    if not os.path.exists(tra_db_dir):
        os.makedirs(tra_db_dir)
        wget_fa_cmd = f'wget https://github.com/DEpt-metagenom/TraDb/raw/main/tra_db/Tra_db_sequences_aa.fa.gz -O {tra_amino_acids}.gz && gunzip {tra_amino_acids}.gz'
        wget_tsv_cmd = f'wget https://github.com/DEpt-metagenom/TraDb/raw/main/tra_db/Tra_db_taxa_positons_aa.tsv.gz -O {tra_amino_data}.gz && gunzip {tra_amino_acids}.gz'
        os.system(wget_fa_cmd)
        os.system(wget_tsv_cmd)

    blastp_cmd = f'blastp -db {amino_acids_db} -query {tra_amino_acids} -num_threads {num_proc} -max_hsps 1 -max_target_seqs 5 -outfmt \'6 qseqid sseqid pident nident qcovhsp evalue gaps slen qlen score\' -subject_besthit | awk  \'$3 >= {mai} && $5 >= {mac} && $6 <= {mae}\' | sed \'1i cluster\thit\tpercent_ident\tnum_ident\tquery_coverage\te-value\tgaps\thit_length\treference_length\tscore\' > {outdir}/{basename}_aa_raw_hits.tsv'

    os.system(blastp_cmd)

    raw_aa_hits = f'{outdir}/{basename}_aa_raw_hits.tsv'
    filtered_aa_hits = f'{outdir}/{basename}_aa_filtered_hits.tsv'
    hit_aa_data = f'{outdir}/{basename}_aa_hit_data.tsv'

    with open(raw_aa_hits, 'r') as fp:
        bp_raw_lines = len(fp.readlines())

    if bp_raw_lines == 0:
        print("No BLASTP hits could be found.")
        exit()

    filter_hits(raw_aa_hits, filtered_aa_hits)

    with open(filtered_aa_hits, 'r') as fp:
        bp_filt_lines = len(fp.readlines())

    if bp_filt_lines == 0:
        print("Filtering removed all BLASTP hits. Check the table of raw hits and try running again with updated criteria.")
        exit()

    get_metadata(filtered_aa_hits, tra_amino_data, hit_aa_data)

    orf_aa_fasta = f'{outdir}/{basename}_aa.fa'
    out_aa_fasta = f'{outdir}/{basename}_aa_hit_sequences.fa'

    filter_fasta(orf_aa_fasta, hit_aa_data, out_aa_fasta)

def run_blastn(dna_db, tra_db_dir, tra_dna, tra_dna_data, basename, outdir, num_proc, mdi, mdc, mde):
    """function to run blastn"""

    if not os.path.exists(f'{dna_db}.ndb'):
        make_dna_db_cmd = f'makeblastdb -in {outdir}/{basename}_dna.fa -dbtype nucl -out {dna_db} 1> /dev/null'
        os.system(make_dna_db_cmd)
        print("Database of ORFs in nucleotide format created")

    if not os.path.exists(tra_db_dir):
        os.makedirs(tra_db_dir)
        wget_fa_cmd = f'wget https://github.com/DEpt-metagenom/TraDb/raw/main/tra_db/Tra_db_sequences_dna.fa.gz -O {tra_dna}.gz && gunzip {tra_dna}.gz'
        wget_tsv_cmd = f'wget https://github.com/DEpt-metagenom/TraDb/raw/main/tra_db/Tra_db_taxa_positons_dna.tsv.gz -O {tra_dna_data}.gz && gunzip {tra_dna}.gz'
        os.system(wget_fa_cmd)
        os.system(wget_tsv_cmd)

    blastn_cmd = f'blastn -db {dna_db} -query {tra_dna} -num_threads {num_proc} -max_hsps 1 -max_target_seqs 5 -outfmt \'6 qseqid sseqid pident nident qcovhsp evalue gaps slen qlen score\' -subject_besthit | awk  \'$3 >= {mdi} && $5 >= {mdc} && $6 <= {mde}\' | sed \'1i cluster\thit\tpercent_ident\tnum_ident\tquery_coverage\te-value\tgaps\thit_length\treference_length\tscore\' > {outdir}/{basename}_dna_raw_hits.tsv'
    os.system(blastn_cmd)

    raw_dna_hits = f'{outdir}/{basename}_dna_raw_hits.tsv'
    filtered_dna_hits = f'{outdir}/{basename}_dna_filtered_hits.tsv'
    hit_dna_data = f'{outdir}/{basename}_dna_hit_data.tsv'

    with open(raw_dna_hits, 'r') as fp:
        bn_raw_lines = len(fp.readlines())

    if bn_raw_lines == 0:
        print("No BLASTP hits could be found.")
        exit()

    filter_hits(raw_dna_hits, filtered_dna_hits)

    with open(filtered_dna_hits, 'r') as fp:
        bn_filt_lines = len(fp.readlines())

    if bn_filt_lines == 0:
        print("Filtering removed all BLASTN hits. Check the table of raw hits and try running again with updated criteria.")
        exit()

    get_metadata(filtered_dna_hits, tra_dna_data, hit_dna_data)

    orf_dna_fasta = f'{outdir}/{basename}_dna.fa'
    out_dna_fasta = f'{outdir}/{basename}_dna_hit_sequences.fa'

    filter_fasta(orf_dna_fasta, hit_dna_data, out_dna_fasta)

def main():
    """main function, parses CLI arguments and uses a basic logic to run the pipeline"""
    print("""
                                          _____          ____  _
           ___  ___ _ __ ___  ___ _ __   |_   _| __ __ _|  _ \| |__
          / __|/ __| '__/ _ \/ _ \ '_ \    | || '__/ _` | | | | '_ \\
          \__ \ (__| | |  __/  __/ | | |   | || | | (_| | |_| | |_) |
          |___/\___|_|  \___|\___|_| |_|   |_||_|  \__,_|____/|_.__/ v0.1
    """)

    # parse arguments
    parser = argparse.ArgumentParser(description="Screen conjugational transfer genes in plasmid sequences.")
    parser.add_argument("input_fasta", help="Input fasta file. Mandatory input.")
    parser.add_argument("--type", "-t", default='both', choices=['dna', 'aa', 'both'], help="Type of search to be run.")
    parser.add_argument("--reading_frames", "-ORF", default='False', action='store_true', help="If the input consits of ORFs. Use this only if you alredy predicted the ORFs of the input. If set, you should also consider to change --type appropriately. Default is False.")
    parser.add_argument("--plasmid_prediction", "-p", default='False', action='store_true', help="If plasmid prediction of the input should be attempted. The output of platon containing the plasmids will be screened for transfer genes")
    parser.add_argument("--number_of_processors", "-np", type=int, default='1', help="Number of processes invoked for the search.")
    parser.add_argument("--output_directory", "-o", default="tra_out", help="Output directory to store results. Output files will be named using the basename of the input file.")
    parser.add_argument("--platon_db_dir", "-pdb", default=f"{home_dir}/platon_db", help="Where the reference database of platon can be found. Default is $HOME/platon_db. If it can not be found automatic install will be attempted.")
    parser.add_argument("--transfer_gene_db_dir", "-tdb", default=f"{home_dir}/tra_db", help="Where the reference conjugational transfer gene database can be found. Default is $HOME/tra_db. If it can not be found automatic install will be attempted.")
    parser.add_argument("--min_aa_identity", "-mai", default='80', help="Minimum percent of identity when screening amino acids (0-100). Default is 80.")
    parser.add_argument("--min_aa_coverage", "-mac", default='80', help="Minimum percent of identity when screening amino acids (0-100). Default is 80.")
    parser.add_argument("--min_aa_evalue", "-mae", default='1e-50', help="Minimum e-value when screening amino acids. Default is 1e-50.")
    parser.add_argument("--min_dna_identity", "-mdi", default='80', help="Minimum percent of identity when screening DNA sequences (0-100). Default is 80.")
    parser.add_argument("--min_dna_coverage", "-mdc", default='80', help="Minimum percent of identity when screening DNA sequences (0-100). Default is 80.")
    parser.add_argument("--min_dna_evalue", "-mde", default='1e-50', help="Minimum e-value when screening DNA sequences. Default is 1e-50.")
    parser.add_argument("--version", "-V", action='version', version=version_number, help="Output version number and exit")
    args = parser.parse_args()

    infile = args.input_fasta
    run_type = args.type
    basename = os.path.splitext(os.path.basename(infile))[0]
    outdir = args.output_directory
    input_orf = args.reading_frames
    plasmid_predict = args.plasmid_prediction
    num_proc = args.number_of_processors
    platon_db_dir = args.platon_db_dir
    tra_db_dir = args.transfer_gene_db_dir
    amino_acids_db = f'{outdir}/{basename}_aa_db'
    tra_amino_acids = f'{tra_db_dir}/Tra_db_sequences_aa.fa'
    tra_amino_data = f'{tra_db_dir}/Tra_db_taxa_positons_aa.tsv'
    mai = args.min_aa_identity
    mac = args.min_aa_coverage
    mae = args.min_aa_evalue
    dna_db = f'{outdir}/{basename}_dna_db'
    tra_dna = f'{tra_db_dir}/Tra_db_sequences_dna.fa'
    tra_dna_data = f'{tra_db_dir}/Tra_db_taxa_positons_dna.tsv'
    mdi = args.min_dna_identity
    mdc = args.min_dna_coverage
    mde = args.min_dna_evalue

    # TODO improve logging
    print(f'Reading input file: {args.input_fasta}')
    print(f'Type of search is set to: {args.type}')
    #print(args.reading_frames)
    #print(args.plasmid_prediction)
    print(f'The run will use {args.number_of_processors} threads')
    print(f'Output will be saved to {args.output_directory}')
    #print(args.platon_db_dir)

    dependencies = ['prodigal', 'bedtools', 'platon', 'wget', 'makeblastdb', 'blastp', 'blastn']

    for dep in dependencies:
        is_tool = os.system(f'which {dep}')
        if is_tool > 0:
            print(f'{dep} is not available')
            exit()

    if input_orf != "False":
        if plasmid_predict == input_orf:
            print("--reading_frames and --plasmid_prediction are both set \"True\". Check your input and run screen_tradb again.")
            exit()

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # logic
    if input_orf == "False" and plasmid_predict == "False":
        run_prodigal(infile, basename, outdir)
        if run_type == "both":
            run_blastn(dna_db, tra_db_dir, tra_dna, tra_dna_data, basename, outdir, num_proc, mdi, mdc, mde)
            run_blastp(amino_acids_db, tra_db_dir, tra_amino_acids, tra_amino_data, basename, outdir, num_proc, mai, mac, mae)
        if run_type == "dna":
            run_blastn(dna_db, tra_db_dir, tra_dna, tra_dna_data, basename, outdir, num_proc, mdi, mdc, mde)
        if run_type == "aa":
            run_blastp(amino_acids_db, tra_db_dir, tra_amino_acids, tra_amino_data, basename, outdir, num_proc, mai, mac, mae)

    elif input_orf != "False" and plasmid_predict == "False":
        if run_type == "dna":
            cp_cmd = f'cp {infile} {outdir}/{basename}_dna.fa'
            os.system(cp_cmd)
            run_blastn(dna_db, tra_db_dir, tra_dna, tra_dna_data, basename, outdir, num_proc, mdi, mdc, mde)
        if run_type == "aa":
            cp_cmd = f'cp {infile} {outdir}/{basename}_aa.fa'
            os.system(cp_cmd)
            run_blastp(amino_acids_db, tra_db_dir, tra_amino_acids, tra_amino_data, basename, outdir, num_proc, mai, mac, mae)

    elif input_orf == "False" and plasmid_predict != "False":
        run_platon(infile, basename, outdir, num_proc, platon_db_dir)
        prodigal_infile = f'{basename}_platon/{basename}.plasmid.fasta'
        run_prodigal(prodigal_infile, basename, outdir)
        if run_type == "both":
            run_blastn(dna_db, tra_db_dir, tra_dna, tra_dna_data, basename, outdir, num_proc, mdi, mdc, mde)
            run_blastp(amino_acids_db, tra_db_dir, tra_amino_acids, tra_amino_data, basename, outdir, num_proc, mai, mac, mae)
        if run_type == "dna":
            run_blastn(dna_db, tra_db_dir, tra_dna, tra_dna_data, basename, outdir, num_proc, mdi, mdc, mde)
        if run_type == "aa":
            run_blastp(amino_acids_db, tra_db_dir, tra_amino_acids, tra_amino_data, basename, outdir, num_proc, mai, mac, mae)

if __name__ == "__main__":
    main()
