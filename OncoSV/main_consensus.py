#!/usr/bin/env python3
import numpy as np
import pandas as pd
from .process_vcf_to_dataframe import process_vcf_to_dataframe
from .consensus_calling import consensus_calling
from .filter_consensus_calls import filter_consensus_calls
from .shared_reads_sv import process_shared_reads
from .header_combine import combine_vcf_lines
from .prepare_vcf_output_file import generate_vcf_from_dataframe

def run_consensus(args):
    if args.chrom == 'all':
        chroms = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')
    
    apply_af_filtering = True  # Default value
    if args.apply_af_filtering is not None:
        apply_af_filtering = args.apply_af_filtering.lower() == "true"

    print("Processing Sniffles VCF file...")
    sniffles_df = process_vcf_to_dataframe(
        args.sniffles,
        chroms,
        qual=args.quality_threshold,
        vcf_format='sniffles',
        lower_sv_size=args.minimum_sv_size,
        upper_sv_size=args.maximum_sv_size,
        sample_id=getattr(args, 'sample_id', None),
        apply_af_filtering=apply_af_filtering
    )
    print(f"Number of variants in Sniffles VCF file: {len(sniffles_df)}")

    print("Processing CuteSV VCF file...")
    cutesv_df = process_vcf_to_dataframe(
        args.cutesv,
        chroms,
        qual=args.quality_threshold,
        vcf_format='cutesv',
        lower_sv_size=args.minimum_sv_size,
        upper_sv_size=args.maximum_sv_size,
        sample_id=getattr(args, 'sample_id', None),
        apply_af_filtering=apply_af_filtering
    )
    print(f"Number of variants in CuteSV VCF file: {len(cutesv_df)}")

    print("Processing SVIM VCF file...")
    svim_df = process_vcf_to_dataframe(
        args.svim,
        chroms,
        qual=args.quality_threshold,
        vcf_format='svim',
        lower_sv_size=args.minimum_sv_size,
        upper_sv_size=args.maximum_sv_size,
        sample_id=getattr(args, 'sample_id', None),
        apply_af_filtering=apply_af_filtering
    )
    print(f"Number of variants in SVIM VCF file: {len(svim_df)}")

    print("Generating consensus calls...")
    consensus_df = consensus_calling(sniffles_df, cutesv_df, svim_df, chroms=chroms)
    consensus_filtered = filter_consensus_calls(consensus_df)
    
    print(f"Number of variants in Consensus VCF file: {len(consensus_filtered)}")

    print("Combining header lines...")
    combined_contigs = combine_vcf_lines(
        args.sniffles,
        args.cutesv,
        '##contig=<ID=',
        chroms,
        extended_chroms=True
    )
    combined_filters = combine_vcf_lines(
        args.sniffles,
        args.cutesv,
        '##FILTER=<ID='
    )

    output_filename = args.out_file
    if args.compress:
        if not output_filename.endswith('.vcf.gz'):
            output_filename += '.vcf.gz' if not output_filename.endswith('.vcf') else '.gz'
        is_compressed = True
    else:
        output_filename += '.vcf' if not output_filename.endswith('.vcf') else ''
        is_compressed = False

    print("Generating VCF output...")
    generate_vcf_from_dataframe(
        consensus_filtered,
        combined_contigs,
        combined_filters,
        output_filename,
        is_compressed
    )

    print(f"VCF file written to {output_filename}")
    print("Consensus structural variant calling completed")
