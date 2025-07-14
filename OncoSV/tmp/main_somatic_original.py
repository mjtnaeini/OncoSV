#!/usr/bin/env python3

import os
from .process_vcf_to_dataframe import process_vcf_to_dataframe
from .identify_variants_withID_proximity import identify_variants
from .prepare_vcf_output_file import generate_vcf_variants

def run_pair(args):
    if args.chrom == 'all':
        chroms = ['chr' + str(i + 1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')

    print("Processing VCF files...")

    tumour_df = process_vcf_to_dataframe(
        args.tumour_consensus,
        chroms,
        qual=args.quality_threshold,
        vcf_format=args.vcf_format,
        lower_sv_size=args.minimum_sv_size,
        upper_sv_size=args.maximum_sv_size,
        sample_id=args.tumour_id
    )

    normal_df = process_vcf_to_dataframe(
        args.normal_consensus,
        chroms,
        qual=0,
        vcf_format=args.vcf_format,
        lower_sv_size=args.minimum_sv_size,
        upper_sv_size=args.maximum_sv_size,
        sample_id=args.normal_id,
        apply_af_filtering=False
    )

    print(f"Identifying somatic and germline variants for {args.svcaller} outputs ...")
    somatic_tumour_df, germline_tumour_df, germline_normal_df, other_normal_df = identify_variants(
        tumour_df,
        normal_df,
        chroms,
        window_size=20
    )

    print(f"Number of somatic structural variants: {len(somatic_tumour_df)}")
    print(f"Number of germline structural variants: {len(germline_tumour_df)}")
    print(f"Number of evidence of germline structural variants in normal sample: {len(germline_normal_df)}")
    print(f"Number of mosaic-normal structural variants in normal sample: {len(other_normal_df)}")

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

