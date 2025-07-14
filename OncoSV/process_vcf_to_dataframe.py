#!/usr/bin/env python3

import pysam
import pandas as pd
import re

def process_sample_data(sample, vcf_format):
    genotype = sample['GT']

    if vcf_format in ['sniffles', 'cutesv', 'consensus']:
        reference_reads = sample.get('DR', 0)
        variant_reads = sample.get('DV', 0)
        genotype_quality = sample.get('GQ', 0)

        reference_reads = int(reference_reads) if reference_reads is not None else 0
        variant_reads = int(variant_reads) if variant_reads is not None else 0

        # Calculate AF for cutesv and sniffles if they use similar structure for DR and DV
        af = variant_reads / (variant_reads + reference_reads) if (variant_reads + reference_reads) > 0 else 0

    elif vcf_format == 'svim':
        svtype = sample.get('SVTYPE', '.')
        # Try to get DP (Total Read Depth), default to 0 if missing
        read_depth = int(sample.get('DP', 0) or 0)

        # Handle special case for BND variants (Breakends)
        if svtype == 'BND':
            reference_reads = 0  # No reference reads for BND
            variant_reads = int(sample.get('SUPPORT', 0) or 0)  # Use SUPPORT
            af = '.'  # No meaningful AF since ref reads are missing

        else:
            # Extract AD (Allelic Depth) if available
            ad_values = sample.get('AD')

            if ad_values is None or not isinstance(ad_values, (list, tuple)) or len(ad_values) != 2:
                # If AD is missing, use SUPPORT for variant_reads
                reference_reads = 0  # No AD, assume no reference reads
                variant_reads = int(sample.get('SUPPORT', 0) or 0)  # Use SUPPORT
                af = '.'  # AF is missing
            else:
                # AD exists, extract reference & variant reads
                reference_reads = int(ad_values[0]) if ad_values[0] is not None else 0
                variant_reads = int(ad_values[1]) if ad_values[1] is not None else 0

                # Compute AF only if reference_reads is available
                af = variant_reads / (reference_reads + variant_reads) if (reference_reads + variant_reads) > 0 else '.'

        genotype_quality = sample.get('GQ', '.')  # Default GQ


    else:
        raise ValueError(f"Unsupported vcf_format: {vcf_format}")  # Catch unexpected cases

    call_record = {
        'Genotype': genotype,
        'GenotypeQuality': genotype_quality,
        'ReferenceReads': reference_reads,
        'VariantReads': variant_reads,
        'AF': af  # Apply AF calculation to all formats where applicable
    }

    return call_record

def convert_svlen(value):
    if isinstance(value, (list, tuple)) and len(value) > 0:
        return int(value[0])
    else:
        try:
            return int(value)
        except:
            return 0

def process_vcf_to_dataframe(vcf_file, chromosomes, qual, vcf_format, lower_sv_size=50, upper_sv_size=1000000, sample_id=None, apply_af_filtering=True):
    processed_data = []

    # Open VCF or compressed VCF file
    with pysam.VariantFile(vcf_file, 'r') as vcf_reader:
        if not sample_id:
            sample_id = list(vcf_reader.header.samples)[0] if vcf_reader.header.samples else 'DefaultSample'

        for record in vcf_reader:
            CHROM = record.chrom

            if CHROM in chromosomes:
                POS = record.pos
                ID = record.id
                REF = record.ref
                ALT = record.alts
                QUAL = record.qual
                FILTER = record.filter.keys()[0] if record.filter else '.'

                info_dict = dict(record.info.items())

                if vcf_format in ['sniffles', 'cutesv', 'consensus']:
                    if info_dict.get('PRECISE', False):
                        filter_status = 'PRECISE'
                    elif info_dict.get('IMPRECISE', False):
                        filter_status = 'IMPRECISE'
                    else:
                        filter_status = '.'
                else:
                    filter_status = '.'

                info_record = {
                    'CHROM': CHROM,
                    'POS': POS,
                    'ID': ID,
                    'REF': REF,
                    'ALT': ALT,
                    'QUAL': QUAL,
                    'FILTER': FILTER,
                    'TYPE': filter_status,
                    'END': record.stop
                }

                info_record.update(info_dict)

                if 'SVTYPE' in info_dict and info_dict['SVTYPE'] == 'BND':
                    if ALT:
                        alt_str = ALT[0]
                        match = re.search(r'N?\[?(chr[\dXY]+):(\d+)\]?N?', alt_str)
                        if match:
                            chr2 = match.group(1)
                            end = int(match.group(2))
                        else:
                            chr2 = '.'
                            end = '.'
                        info_record['CHR2'] = chr2
                        info_record['END'] = end
                    else:
                        chr2 = '.'
                        end = '.'
                        info_record['CHR2'] = chr2
                        info_record['END'] = end

                sample_data = record.samples.get(sample_id)
                if sample_data is not None:
                    call_record = process_sample_data(sample_data, vcf_format)
                    final_record = {**info_record, **call_record, 'Sample': sample_id}
                    processed_data.append(final_record)

    processed_df = pd.DataFrame(processed_data)

    column_order = list(processed_df.columns)
    column_order.remove('Sample')
    column_order.append('Sample')
    processed_df = processed_df[column_order]

    processed_df['CHR2'] = processed_df['CHR2'].fillna(processed_df['CHROM'])
    processed_df.rename(columns={'CHR2': 'CHROM2'}, inplace=True)

    processed_df['END'] = processed_df['END'].fillna(processed_df['POS'])
    processed_df['SVLEN'] = processed_df['SVLEN'].apply(convert_svlen)

    filtered_df = processed_df[
    (processed_df['SVTYPE'].isin(['BND', 'INV'])) |
    ((processed_df['SVLEN'].abs() >= lower_sv_size) & (processed_df['SVLEN'].abs() <= upper_sv_size))
]
    filtered_df = filtered_df[filtered_df['QUAL'] >= qual]

    # Apply AF filtering only if enabled
    if apply_af_filtering:
        # Separate BND and non-BND variants
        bnd_variants = filtered_df[filtered_df['SVTYPE'] == 'BND']
        non_bnd_variants = filtered_df[filtered_df['SVTYPE'] != 'BND']

        # Apply AF filtering only to non-BND variants
        non_bnd_variants = non_bnd_variants.dropna(subset=['AF'])
        non_bnd_variants = non_bnd_variants.loc[
            non_bnd_variants.groupby(['CHROM', 'POS', 'END'])['AF'].idxmax()
        ]
        # Concatenate BND and filtered non-BND variants
        filtered_df = pd.concat([bnd_variants, non_bnd_variants])

    return filtered_df

