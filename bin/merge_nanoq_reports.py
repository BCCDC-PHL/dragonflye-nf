#!/usr/bin/env python3

import argparse
import collections
import csv
import json
import sys

def parse_nanoq_report(nanoq_report_path):
    fieldname_lookup = {
        'reads':          'total_reads',
        'bases':          'total_bases',
        'n50':            'read_n50',
        'shortest':       'shortest_read_length',
        'longest':        'longest_read_length',
        'mean_length':    'mean_read_length',
        'median_length':  'median_read_length',
        'mean_quality':   'mean_base_quality',
        'median_quality': 'median_base_quality',
    }
    nanoq_report = []
    with open(nanoq_report_path, 'r') as f:
        reader = csv.DictReader(f, dialect='unix')
        for row in reader:
            nanoq_record = {}
            for field in row:
                if field in fieldname_lookup:
                    fieldname = fieldname_lookup[field]
                    nanoq_record[fieldname] = row[field]
                else:
                    nanoq_record[field] = row[field]
                    
            nanoq_report.append(nanoq_record)

    return nanoq_report


def main(args):
    nanoq_pre_filter = { k + '_before_filtering': v for (k, v) in parse_nanoq_report(args.pre_filter)[0].items() }
    nanoq_post_filter = { k + '_after_filtering': v for (k, v) in parse_nanoq_report(args.post_filter)[0].items() }

    nanoq_merged = nanoq_pre_filter.copy()
    nanoq_merged.update(nanoq_post_filter)
    
    if args.sample_id:
        output_fieldnames = ['sample_id']
        nanoq_merged['sample_id'] = args.sample_id
    else:
        output_fieldnames = []

    output_fieldnames += [
        "total_reads_before_filtering",
        "total_reads_after_filtering",
        "total_bases_before_filtering",
        "total_bases_after_filtering",
        "mean_read_length_before_filtering",
        "mean_read_length_after_filtering",
        "median_read_length_before_filtering",
        "median_read_length_after_filtering",
        "shortest_read_length_before_filtering",
        "shortest_read_length_after_filtering",
        "longest_read_length_before_filtering",
        "longest_read_length_after_filtering",
        "read_n50_before_filtering",
        "read_n50_after_filtering",
        "mean_base_quality_before_filtering",
        "mean_base_quality_after_filtering",
        "median_base_quality_before_filtering",
        "median_base_quality_after_filtering",
    ]

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    writer.writerow(nanoq_merged)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-id')
    parser.add_argument('--pre-filter')
    parser.add_argument('--post-filter')
    args = parser.parse_args()
    main(args)
    
