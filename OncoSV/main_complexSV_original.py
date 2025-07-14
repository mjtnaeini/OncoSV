#!/usr/bin/env python3

import pandas as pd
import os
import argparse
from .process_vcf_to_dataframe import process_vcf_to_dataframe
from .shared_reads_sv import process_shared_reads, process_sv_data_with_sv_count, process_breakpoints, add_overlapping_column
from .find_network_sv import group_by_group, identify_networks

def run_complexSV(args):
    if args.chrom == 'all':
        chroms = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')

    print("Processing Input VCF file...")
    vcf = process_vcf_to_dataframe(
            args.vcf,
            chroms,
            qual=args.qual,
            vcf_format=args.vcf_format,
            lower_sv_size=args.minimum_sv_size,
            upper_sv_size=args.maximum_sv_size,
            sample_id=getattr(args, 'sample_id', None),
            apply_af_filtering=False)
    print(f"Number of variants processed: {len(vcf)}")

    print("Processing shared reads...")
    shared_reads = process_shared_reads(vcf)
    print(f"Number of shared reads identified: {len(shared_reads)}")

    print("Processing SV shared reads with counts...")
    shared_sv_counts = process_sv_data_with_sv_count(shared_reads)
    shared_sv_counts_breakopints = process_breakpoints(shared_sv_counts)
    shared_sv_counts_breakopints_overlap = add_overlapping_column(shared_sv_counts_breakopints)
    shared_sv_counts_breakopints_overlap_path = os.path.join(args.output_dir, f"{args.label_prefix + '_' if args.label_prefix else ''}shared_sv_counts_breakopints_overlap.csv")
    print(f"Saving SV shared reads with counts and overlapping breakpoints to {shared_sv_counts_breakopints_overlap_path}...")
    print(f"Number of shared SV counts: {len(shared_sv_counts_breakopints_overlap)}")
    shared_sv_counts_breakopints_overlap.to_csv(shared_sv_counts_breakopints_overlap_path, index=False)

    print("Grouping by complex sv groups...")
    complexSV_df = group_by_group(shared_sv_counts_breakopints_overlap)
    
    # Find number of complex SV groups
    all_groups = []

    for group in complexSV_df['group'].astype(str):
        split_groups = group.split(';')
        for item in split_groups:
            all_groups.append(int(item))

    # Find the maximum number
    max_group_number = max(all_groups)
    print(f"Number of complex SV groups: {max_group_number}")

    print("Identifying networks...")
    complexSV_network_df = identify_networks(complexSV_df)
    complexSV_networks_path = os.path.join(args.output_dir, f"{args.label_prefix + '_' if args.label_prefix else ''}complexSV_groups_networks.csv")
    print(f"Saving complex SV group and networks data to {complexSV_networks_path}...")
    unique_network_count = complexSV_network_df['Network'].nunique()
    print(f"Number of complex SV networks: {unique_network_count}")
    complexSV_network_df.to_csv(complexSV_networks_path, index=False)

    print("Complex SV network calling completed successfully.")

