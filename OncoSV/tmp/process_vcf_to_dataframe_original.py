#!/usr/bin/env python3

import pysam
import pandas as pd
import re

def process_sample_data(sample, vcf_format):
    genotype = sample['GT']

    if vcf_format in ['sniffles', 'cutesv', 'consensus']:
        reference_reads = sample['DR']
        variant_reads = sample['DV']
        genotype_quality = sample['GQ']
        # Calculate AF for cutesv and sniffles if they use similar structure for DR and DV
        af = variant_reads / (variant_reads + reference_reads) if (variant_reads + reference_reads) > 0 else 0

    elif vcf_format == 'svim':
        ad_values = sample.get('AD', (0, 0))  # Default to (0, 0) if missing or None

        if isinstance(ad_values, (list, tuple)):
            if len(ad_values) == 2:
                # Safely handle None values in AD
                reference_reads = int(ad_values[0]) if ad_values[0] is not None else 0
                variant_reads = int(ad_values[1]) if ad_values[1] is not None else 0
            else:
                reference_reads, variant_reads = 0, 0  # Default to 0 if AD length is invalid
        else:
            reference_reads, variant_reads = 0, 0  # Default to 0 if AD is not a list/tuple

        genotype_quality = sample.get('GQ', '.')
        af = variant_reads / (variant_reads + reference_reads) if (variant_reads + reference_reads) > 0 else 0

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

def process_vcf_to_dataframe(vcf_file, chromosomes, qual, vcf_format, lower_sv_size=50, upper_sv_size=1000000, sample_id=None):
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

    # Retain only the SVs with the highest AF for exact coordinates
    filtered_df['AF'] = pd.to_numeric(filtered_df['AF'], errors='coerce')
    filtered_df = filtered_df.dropna(subset=['AF'])  # Remove rows with NaN AF values
    filtered_df = filtered_df.loc[
        filtered_df.groupby(['CHROM', 'POS', 'END'])['AF'].idxmax()
    ]

    return filtered_df

