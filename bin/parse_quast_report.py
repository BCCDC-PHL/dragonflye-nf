#!/usr/bin/env python3

import argparse
import collections
import csv
import json
import sys


def parse_transposed_quast_report(transposed_quast_report_path):
    """
    """
    field_lookup = collections.OrderedDict()
    field_lookup['Assembly'] = 'assembly_id'
    field_lookup['Total length'] = 'total_length'
    field_lookup['# contigs'] = 'num_contigs'
    field_lookup['Largest contig'] = 'largest_contig'
    field_lookup['N50'] = 'assembly_N50'
    field_lookup['N75'] = 'assembly_N75'
    field_lookup['L50'] = 'assembly_L50'
    field_lookup['L75'] = 'assembly_L75'
    field_lookup["# N's per 100 kbp"] = 'num_N_per_100_kb'
    field_lookup['# contigs (>= 0 bp)'] = 'num_contigs_gt_0_bp'
    field_lookup['# contigs (>= 1000 bp)'] = 'num_contigs_gt_1000_bp'
    field_lookup['# contigs (>= 5000 bp)'] = 'num_contigs_gt_5000_bp'
    field_lookup['# contigs (>= 10000 bp)'] = 'num_contigs_gt_10000_bp'
    field_lookup['# contigs (>= 25000 bp)'] = 'num_contigs_gt_25000_bp'
    field_lookup['# contigs (>= 50000 bp)'] = 'num_contigs_gt_50000_bp'
    field_lookup['Total length (>= 0 bp)'] = 'total_length_gt_0_bp'
    field_lookup['Total length (>= 1000 bp)'] = 'total_length_gt_1000_bp'
    field_lookup['Total length (>= 5000 bp)'] = 'total_length_gt_5000_bp'
    field_lookup['Total length (>= 10000 bp)'] = 'total_length_gt_10000_bp'
    field_lookup['Total length (>= 25000 bp)'] = 'total_length_gt_25000_bp'
    field_lookup['Total length (>= 50000 bp)'] = 'total_length_gt_50000_bp'


    int_fields = [
        'total_length',
        'num_contigs',
        'largest_contig',
        'assembly_N50',
        'assembly_N75',
        'assembly_L50',
        'assembly_L75',
        'num_contigs_gt_0_bp',
        'num_contigs_gt_1000_bp',
        'num_contigs_gt_5000_bp',
        'num_contigs_gt_10000_bp',
        'num_contigs_gt_25000_bp',
        'num_contigs_gt_50000_bp',
        'total_length_gt_0_bp',
        'total_length_gt_1000_bp',
        'total_length_gt_5000_bp',
        'total_length_gt_10000_bp',
        'total_length_gt_25000_bp',
        'total_length_gt_50000_bp',
    ]

    float_fields = [
        'num_N_per_100_kb',
    ]

    parsed_report = []
    with open(transposed_quast_report_path, 'r', newline='') as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        for row in reader:
            r = collections.OrderedDict()
            for f in field_lookup:
                r[field_lookup[f]] = row[f]

            for f in int_fields:
                try:
                    r[f] = int(r[f])
                except ValueError as e:
                    r[f] = None

            for f in float_fields:
                try:
                    r[f] = float(r[f])
                except ValueError as e:
                    r[f] = None

            parsed_report.append(r)

    return parsed_report



def main():

    
    parser = argparse.ArgumentParser()
    parser.add_argument('transposed_quast_report')
    args = parser.parse_args()

    output_fieldnames = [
        'assembly_id',
        'total_length',
        'num_contigs',
        'largest_contig',
        'assembly_N50',
        'assembly_N75',
        'assembly_L50',
        'assembly_L75',
        'num_contigs_gt_0_bp',
        'num_contigs_gt_1000_bp',
        'num_contigs_gt_5000_bp',
        'num_contigs_gt_10000_bp',
        'num_contigs_gt_25000_bp',
        'num_contigs_gt_50000_bp',
        'total_length_gt_0_bp',
        'total_length_gt_1000_bp',
        'total_length_gt_5000_bp',
        'total_length_gt_10000_bp',
        'total_length_gt_25000_bp',
        'total_length_gt_50000_bp',
        'num_N_per_100_kb',
    ]

    report = parse_transposed_quast_report(args.transposed_quast_report)
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames)
    writer.writeheader()
    for record in report:
        writer.writerow(record)


if __name__ == '__main__':
    main()
