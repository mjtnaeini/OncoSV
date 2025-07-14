#!/usr/bin/env python3

import os
import pandas as pd
from .process_vcf_to_dataframe import process_vcf_to_dataframe
from .identify_variants_withID_proximity import identify_variants
from .prepare_vcf_output_file import generate_vcf_variants

def detect_vcf_format(filename):
    """Automatically detects VCF format based on filename."""
    filename = filename.lower()
    if "sniffles" in filename:
        return "sniffles"
    elif "cutesv" in filename:
        return "cutesv"
    elif "svim" in filename:
        return "svim"
    else:
        raise ValueError(f"Unknown VCF format for file: {filename}")

def run_pair(args):
    if args.chrom == 'all':
        chroms = ['chr' + str(i + 1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')

    print("Processing tumour VCF file...")
    tumour_df = process_vcf_to_dataframe(
        args.tumour_consensus,
        chroms,
        qual=args.quality_threshold,
        vcf_format=args.vcf_format,
        lower_sv_size=args.minimum_sv_size,
        upper_sv_size=args.maximum_sv_size,
        sample_id=args.tumour_id
    )

    if args.normal_mode == "single":
        print("Processing single normal VCF file...")
        normal_df = process_vcf_to_dataframe(
            args.normal_sample,
            chroms,
            qual=0,
            vcf_format=args.vcf_format,
            lower_sv_size=args.minimum_sv_size,
            upper_sv_size=args.maximum_sv_size,
            sample_id=args.normal_id,
            apply_af_filtering=False
        )
    elif args.normal_mode == "multi":
        print("Processing multiple normal VCF files...")
        normal_dfs = []
        for normal_vcf in [args.normal_sample1, args.normal_sample2, args.normal_sample3]:
            vcf_format = detect_vcf_format(normal_vcf)
            df = process_vcf_to_dataframe(
                normal_vcf,
                chroms,
                qual=0,
                vcf_format=vcf_format,
                lower_sv_size=args.minimum_sv_size,
                upper_sv_size=args.maximum_sv_size,
                sample_id=args.normal_id,
                apply_af_filtering=False
            )
            normal_dfs.append(df)
        normal_df = pd.concat(normal_dfs, ignore_index=True)
        print(f"Total variants after merging normal samples: {len(normal_df)}")
        
        # Save merged normal data as CSV if specified as 'true'
        if args.save_merged_normal.lower() == "true":
            merged_csv_path = os.path.join(args.out_dir, "merged_normal_samples.csv")
            normal_df.to_csv(merged_csv_path, index=False)
            print(f"Merged normal samples saved to {merged_csv_path}")
    else:
        raise ValueError("Invalid normal mode. Choose 'single' or 'multi'.")

    print(f"Identifying somatic and germline variants for {args.svcaller} outputs ...")
    somatic_tumour_df, germline_tumour_df, germline_normal_df, other_normal_df = identify_variants(
        tumour_df,
        normal_df,
        chroms,
        window_size=200
    )

    print(f"Number of somatic structural variants: {len(somatic_tumour_df)}")
    print(f"Number of germline structural variants: {len(germline_tumour_df)}")
    print(f"Number of germline variants in normal samples: {len(germline_normal_df)}")
    print(f"Number of mosaic-normal variants: {len(other_normal_df)}")

    # Construct the output filenames
    somatic_output_filename = os.path.join(args.out_dir, f'{args.svcaller}_somatic_variants.vcf')
    germline_tumour_output_filename = os.path.join(args.out_dir, f'{args.svcaller}_germline_variants.vcf')
    germline_normal_evidence_output_filename = os.path.join(args.out_dir, f'{args.svcaller}_germline_normal_evidence.vcf')
    mosaic_normal_output_filename = os.path.join(args.out_dir, f'{args.svcaller}_mosaic_normal.vcf')

    output_directories = [
        os.path.dirname(somatic_output_filename),
        os.path.dirname(germline_tumour_output_filename),
        os.path.dirname(germline_normal_evidence_output_filename),
        os.path.dirname(mosaic_normal_output_filename)
    ]

    for dir in output_directories:
        if not os.path.exists(dir):
            os.makedirs(dir)
            print(f"Created directory: {dir}")

    if args.compress:
        somatic_output_filename += '.gz'
        germline_tumour_output_filename += '.gz'
        germline_normal_evidence_output_filename += '.gz'
        mosaic_normal_output_filename += '.gz'

    if args.only_somatic:
        print("Generating VCF file for somatic tumour variants...")
        generate_vcf_variants(
            somatic_tumour_df,
            args.tumour_consensus,
            somatic_output_filename,
            is_compressed=args.compress,
            svcaller=args.svcaller,
            include_variant_ID=True,
            sample_id=args.patient_id if args.patient_id else None
        )
    else:
        print("Generating VCF files for all variant types...")
        generate_vcf_variants(
            somatic_tumour_df,
            args.tumour_consensus,
            somatic_output_filename,
            is_compressed=args.compress,
            svcaller=args.svcaller,
            include_variant_ID=True,
            sample_id=args.patient_id if args.patient_id else None
        )
        generate_vcf_variants(
            germline_tumour_df,
            args.tumour_consensus,
            germline_tumour_output_filename,
            is_compressed=args.compress,
            svcaller=args.svcaller,
            include_variant_ID=True,
            sample_id=args.patient_id if args.patient_id else None
        )
        generate_vcf_variants(
            germline_normal_df,
            args.normal_consensus,
            germline_normal_evidence_output_filename,
            is_compressed=args.compress,
            svcaller=args.svcaller,
            include_variant_ID=True,
            sample_id=args.patient_id if args.patient_id else None
        )
        generate_vcf_variants(
            other_normal_df,
            args.normal_consensus,
            mosaic_normal_output_filename,
            is_compressed=args.compress,
            svcaller=args.svcaller,
            include_variant_ID=True,
            sample_id=args.patient_id if args.patient_id else None
        )

    print("Somatic variant calling completed")
