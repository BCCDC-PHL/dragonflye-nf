#!/usr/bin/env python3

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys

from pathlib import Path


def reverse_complement(seq: str):
    """
    Given a DNA sequences, this function returns the reverse complement sequence.
    """
    complement_base = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'N': 'N',
        'a': 't',
        'c': 'g',
        'g': 'c',
        't': 'a',
        'n': 'n',
    }

    for idx, base in enumerate(seq):
        if base not in complement_base.keys():
            seq[idx] = 'N'
            
    revcomp_seq = ''.join([complement_base[x] for x in seq][::-1])

    return revcomp_seq


def parse_input_fasta(input_fasta_path: Path):
    """
    Parse fasta file into list of dicts, each with keys: ['id', 'description', 'seq'].
    the value associated with 'description' is a dict with keys corresponding to
    the tags in the fasta description line.
    
    :param input_fasta_path: Path to fasta file.
    :type input_fasta_path: pathlib.Path
    :return: Parsed sequence
    :rtype: list[dict[str, object]]
    """
    parsed_fasta = []
    with open(input_fasta_path, 'r') as f:
        current_contig = {}
        current_description = {}
        current_seq = ""
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq != "":
                    current_contig['seq'] = current_seq
                    parsed_fasta.append(current_contig)
                current_contig = {}
                current_description = {}
                current_seq = ""
                seq_id = line.lstrip('>').split()[0]
                current_contig['id'] = seq_id
                try:
                    description_string = ' '.join(line.lstrip('>').replace(' round(s)', '_round(s)').split()[1:])
                    description_list = description_string.split()
                    for tag in description_list:
                        tag_split = tag.split('=')
                        current_description[tag_split[0]] = tag_split[1]
                    current_contig['description'] = current_description
                except Exception as e:
                    print("error parsing description", e, sys.stderr)
                    pass
            else:
                current_seq += line

        current_contig['seq'] = current_seq
        parsed_fasta.append(current_contig)

    return parsed_fasta


def find_longest_query(query_fasta_path: Path):
    """
    Parse fasta file and find the longest sequence.

    :return: Longest query. Dict with keys: ['id', 'description', 'seq', 'length']
    :rtype: dict
    """
    current_query_length = 0
    longest_query_length = 0
    current_query = {'length': 0}
    longest_query = {'length': 0}
    seq_id = ""
    description = ""

    with open(query_fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if longest_query is None:
                    longest_query = current_query
                elif current_query['length'] >= longest_query['length']:
                    longest_query = current_query
                current_query = {}
                seq_id = line.lstrip('>').split()[0]
                current_query['id'] = seq_id
                try:
                    description = ' '.join(line.lstrip('>').split()[1:])
                    current_query['description'] = description
                except Exception as e:
                    print("Error parsing description for query seq: " + seq_id, sys.stderr)
                    current_query['description'] = ''
            else:
                if 'seq' not in current_query:
                    current_query['seq'] = line
                else:
                    current_query['seq'] += line
                if 'length' not in current_query:
                    current_query['length'] = len(line)
                else:
                    current_query['length'] += len(line)

    return longest_query


def prepare_blast_db(parsed_seq: dict, longest_query_length: int):
    """
    """
    tmp_fasta_filename = parsed_seq['id'].replace(' ', '_') + '.fa'
    db_sequence = parsed_seq['seq']

    # Ensure that we have full contiguous sequence across start/end
    # of the db sequence, sufficient for longest query seq
    dup_length = min(len(db_sequence), longest_query_length)
    db_sequence = db_sequence + db_sequence[:dup_length]

    with open(tmp_fasta_filename, 'w') as f:
        f.write(''.join(['>', parsed_seq['id'], '\n']))
        f.write(db_sequence)
        f.write('\n')

    makeblastdb_cmd = [
        'makeblastdb',
        '-dbtype',
        'nucl',
        '-in',
        tmp_fasta_filename
    ]

    process = subprocess.Popen(makeblastdb_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = process.communicate()
    if err:
        print('makeblastdb encountered an error:\n' + err.decode(), sys.stderr)
        exit(-1)

    return tmp_fasta_filename


def blast_start_genes(blast_db_path: Path, start_genes_fasta: Path, subject_length: int, num_threads=1):
    """
    """

    blast_output = []

    tblastn_cmd = [
        'tblastn',
        '-db', blast_db_path,
        '-query', start_genes_fasta,
        '-outfmt', '6 qseqid sstart send pident qlen qseq qstart qend bitscore',
        '-num_threads', str(num_threads),
    ]

    process = subprocess.Popen(tblastn_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    blast_stdout, blast_stderr = process.communicate()
    process.wait()
    if blast_stderr:
        print('BLAST encountered an error:\n' + blast_err.decode(), sys.stderr)
        exit(-1)

    blast_fields = [
        'query_seq_id',
        'subject_start',
        'subject_end',
        'percent_identity',
        'query_length',
        'query_seq',
        'query_start',
        'query_end',
        'bitscore',
    ]

    int_fields = set([
        'subject_start',
        'subject_end',
        'query_start',
        'query_end',
        'query_length',
    ])

    float_fields = set([
        'percent_identity',
        'bitscore',
    ])

    for line in blast_stdout.decode().splitlines():
        line_split = line.split()
        blast_line = {}
        for idx, field in enumerate(blast_fields):
            value = line_split[idx]
            if field in int_fields:
                try:
                    value = int(value)
                except ValueError as e:
                    value = None
            elif field in float_fields:
                try:
                    value = float(value)
                except ValueError as e:
                    value = None

            if field == 'subject_start' or field == 'query_start' and value is not None:
                value -= 1

            blast_line[field] = value

        blast_line['percent_query_coverage'] = 100.0 * len(blast_line['query_seq']) / blast_line['query_length']
        if blast_line['subject_start'] <= blast_line['subject_end']:
            blast_line['start_pos'] = blast_line['subject_start']
            blast_line['flip'] = False
        else:
            blast_line['start_pos'] = blast_line['subject_start'] + 1
            blast_line['flip'] = True
        if blast_line['start_pos'] >= subject_length:
            blast_line['start_pos'] -= subject_length

        blast_output.append(blast_line)

    return blast_output


def find_best_hit(blast_hits: list, identity_threshold=90.0, coverage_threshold=95.0):
    """
    """

    best_hit = None
    best_bitscore = 0

    for hit in blast_hits:

        if hit['percent_identity'] >= identity_threshold and \
           hit['percent_query_coverage'] >= coverage_threshold and \
           hit['query_start'] == 0 and \
           hit['bitscore'] > best_bitscore:
            best_hit = hit
            best_bitscore = hit['bitscore']
            
    return best_hit


def rotate_sequence(seq: dict, best_hit: dict):
    """
    """
    rotated_description = seq['description']
    rotated_description['rotated'] = "true"
    rotated_description['start_gene_id'] = best_hit['query_seq_id']
    rotated = {
        'id': seq['id'],
        'description': rotated_description,
    }
    start_pos = best_hit['start_pos']
    flip = best_hit['flip']
    unrotated_seq = seq['seq']
    rotated_seq = unrotated_seq[start_pos:] + unrotated_seq[:start_pos]
    rev_comp_rotated_seq = reverse_complement(rotated_seq)

    if flip:
        rotated['seq'] = rev_comp_rotated_seq
    else:
        rotated['seq'] = rotated_seq

    return rotated
    

def main(args):
    """
    """

    if not(os.path.exists(args.start_genes)):
        print("Error: start genes file " + str(args.start_genes) + " does not exist.")
        exit(-1)

    start_genes_abspath = os.path.abspath(args.start_genes)
    tmp_dir_abspath = os.path.abspath(args.tmp)

    longest_query = find_longest_query(start_genes_abspath)

    if 'length' not in longest_query:
        print("Error: failed to find length of longest query seq.", sys.stderr)
        exit()
    longest_query_nucleotide_length = 3 * longest_query['length']

    original_current_working_dir = os.getcwd()
    if not os.path.exists(tmp_dir_abspath):
        os.makedirs(tmp_dir_abspath)

    input_fasta = parse_input_fasta(args.input)

    os.chdir(tmp_dir_abspath)

    rotation_result_fieldnames = [
        'seq_id',
        'seq_length',
        'starting_gene_id',
        'gene_start_position',
        'strand',
        'percent_identity',
        'percent_coverage'
    ]
    rotation_results_writer = csv.DictWriter(sys.stdout, fieldnames=rotation_result_fieldnames,
                                             dialect='unix', quoting=csv.QUOTE_MINIMAL)
    rotation_results_writer.writeheader()

    for contig in input_fasta:
        input_seq_length = len(contig['seq'])
    
        blast_db_path = prepare_blast_db(contig, longest_query_nucleotide_length)
    
        blast_hits = blast_start_genes(blast_db_path, start_genes_abspath, input_seq_length, num_threads=args.blast_threads)

        best_hit = find_best_hit(blast_hits)

        if best_hit is not None:
            rotation_results = {
                'seq_id': contig['id'],
                'seq_length': len(contig['seq']),
                'starting_gene_id': best_hit['query_seq_id'],
                'gene_start_position': best_hit['start_pos'],
                'strand': 'reverse' if best_hit['flip'] else 'forward',
                'percent_identity': '%.1f' % best_hit['percent_identity'],
                'percent_coverage': '%.1f' % best_hit['percent_query_coverage'],
            }
            rotation_results_writer.writerow(rotation_results)

            output_contig = rotate_sequence(contig, best_hit)

        else:
            rotation_results = {
                'seq_id': contig['id'],
                'seq_length': len(contig['seq']),
                'starting_gene_id': '',
                'gene_start_position': '',
                'strand': '',
                'percent_identity': '',
                'percent_coverage': '',
            }
            rotation_results_writer.writerow(rotation_results)

            output_contig = contig

        with open(args.output, 'a') as f:
            description_string = ' '.join([str(x) + '=' + str(y) for x, y in output_contig['description'].items()])
            f.write(''.join(['>', output_contig['id'], ' ', description_string, '\n']))
            for i in range(0, len(output_contig['seq']), args.output_line_length):
                f.write(output_contig['seq'][i:i+args.output_line_length])
                f.write('\n')

    os.chdir(original_current_working_dir)

    src_assembly = Path(tmp_dir_abspath) / args.output.name
    dest_assembly = args.output.parent / args.output.name
    shutil.move(src_assembly, dest_assembly)

    if not args.no_cleanup and os.path.exists(tmp_dir_abspath):        
        shutil.rmtree(tmp_dir_abspath, ignore_errors=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=Path, help="Path to plasmid to be rotated")
    parser.add_argument('-o', '--output', type=Path, default="./rotated.fa", help="Path to write rotated plasmid")
    parser.add_argument('--output-line-length', type=int, default=60, help="Line length for output fasta file.")
    parser.add_argument('-s', '--start-genes', type=Path, help="Path to 'start genes' (dnaA, repA) file (fasta format)")
    parser.add_argument('-t', '--blast-threads', type=int, default=1, help="Number of CPU threads to use for tblastn.")
    parser.add_argument('--tmp', type=Path, default="./tmp", help="Directory to use for temporary files.")
    parser.add_argument('--no-cleanup', action='store_true', help="Do not remove temp directory.")
    args = parser.parse_args()
    main(args)
