#!/usr/bin/env python3

import pandas as pd

def identify_variants(tumour, normal, chromosomes, window_size=200):
    # Save original setting
    original_setting = pd.options.mode.chained_assignment

    # Suppress SettingWithCopyWarning
    pd.options.mode.chained_assignment = None

    # Filter dataframes for specified chromosomes
    tumour = tumour[tumour['CHROM'].isin(chromosomes)].copy()
    normal = normal[normal['CHROM'].isin(chromosomes)].copy()

    # Convert the END and SVLEN columns to numeric, setting errors='coerce' to convert invalid data to NaN
    tumour['END'] = pd.to_numeric(tumour['END'], errors='coerce')
    normal['END'] = pd.to_numeric(normal['END'], errors='coerce')
    tumour['SVLEN'] = pd.to_numeric(tumour['SVLEN'], errors='coerce')
    normal['SVLEN'] = pd.to_numeric(normal['SVLEN'], errors='coerce')

    # Ensure SVLEN is absolute for deletions
    tumour.loc[tumour['SVTYPE'] == 'DEL', 'SVLEN'] = tumour['SVLEN'].abs()
    normal.loc[normal['SVTYPE'] == 'DEL', 'SVLEN'] = normal['SVLEN'].abs()

    # Initialize new columns for variant classification and ID
    tumour['variant_type'] = 'unknown'
    tumour['variant_ID'] = None
    normal['variant_ID'] = None

    # Counter for generating unique numbers
    variant_counter = 1

    # Iterate over each variant in the tumour dataframe
    for index, tumour_row in tumour.iterrows():
        chrom1 = tumour_row['CHROM']
        pos1 = tumour_row['POS']
        chrom2 = tumour_row.get('CHROM2', None)
        pos2 = tumour_row.get('END', None)
        svtype = tumour_row['SVTYPE']
        tumour_svlen = tumour_row['SVLEN']

        # Define the search windows
        window_start1 = pos1 - window_size
        window_end1 = pos1 + window_size
        window_start2 = pos2 - window_size if pos2 else None
        window_end2 = pos2 + window_size if pos2 else None

        # Match logic for different SVTYPEs
        if svtype in ['DUP', 'DEL'] and tumour_svlen > 1000:
            match_normal = normal[
                (normal['CHROM'] == chrom1) & (normal['POS'] >= window_start1) & (normal['POS'] <= window_end1) &
                (normal['SVTYPE'] == svtype)
            ]
            match_reverse = normal[
                (normal['CHROM'] == chrom2) & (normal['POS'] >= window_start2) & (normal['POS'] <= window_end2) &
                (normal['SVTYPE'] == svtype)
            ]
            match = match_normal if not match_normal.empty else match_reverse

        else:
            match_normal = normal[
                (normal['CHROM'] == chrom1) & (normal['POS'] >= window_start1) & (normal['POS'] <= window_end1) &
                (normal['CHROM2'] == chrom2) & (normal['END'] >= window_start2) & (normal['END'] <= window_end2) &
                (normal['SVTYPE'] == svtype)
            ]
            match_reverse = normal[
                (normal['CHROM'] == chrom2) & (normal['POS'] >= window_start2) & (normal['POS'] <= window_end2) &
                (normal['CHROM2'] == chrom1) & (normal['END'] >= window_start1) & (normal['END'] <= window_end1) &
                (normal['SVTYPE'] == svtype)
            ]
            match = match_normal if not match_normal.empty else match_reverse

        # Handle cases with multiple matches
        if not match.empty:
            if len(match) > 1:
                match['proximity'] = match.apply(lambda row: abs(row['POS'] - pos1), axis=1)
                closest_match = match.loc[match['proximity'].idxmin()]
                match = pd.DataFrame([closest_match])

            normal_svlen = match.iloc[0]['SVLEN']
            if svtype == 'INS' and tumour_svlen > 0 and normal_svlen > 0:
                mean_svlen = (tumour_svlen + normal_svlen) / 2
                svlen_sd = ((tumour_svlen - mean_svlen)**2 + (normal_svlen - mean_svlen)**2)**0.5 / 2
                sd_threshold = 0.3 * normal_svlen
                if svlen_sd > sd_threshold:
                    tumour.at[index, 'variant_type'] = 'somatic'
                    tumour.at[index, 'variant_ID'] = f'somatic.{svtype}.{variant_counter}'
                    variant_counter += 1
                    continue

            if svtype in ['DUP', 'DEL'] and tumour_svlen > 1000:
                mean_svlen = (tumour_svlen + normal_svlen) / 2
                svlen_sd = ((tumour_svlen - mean_svlen)**2 + (normal_svlen - mean_svlen)**2)**0.5 / 2
                sd_threshold = 0.2 * normal_svlen
                if svlen_sd > sd_threshold:
                    tumour.at[index, 'variant_type'] = 'somatic'
                    tumour.at[index, 'variant_ID'] = f'somatic.{svtype}.{variant_counter}'
                    variant_counter += 1
                    continue

            match_index = match.index[0]
            existing_id = normal.at[match_index, 'variant_ID']
            if existing_id:
                variant_id = existing_id
            else:
                variant_id = f'germline.{svtype}.{variant_counter}'
                normal.at[match_index, 'variant_ID'] = variant_id
                variant_counter += 1
            tumour.at[index, 'variant_type'] = 'germline-normal'
            tumour.at[index, 'variant_ID'] = variant_id
        else:
            tumour.at[index, 'variant_type'] = 'somatic'
            tumour.at[index, 'variant_ID'] = f'somatic.{svtype}.{variant_counter}'
            variant_counter += 1

    # Split the tumour dataframe into two based on variant_type
    somatic_tumour = tumour[tumour['variant_type'] == 'somatic'].drop(columns=['variant_type'])
    germline_tumour = tumour[tumour['variant_type'] == 'germline-normal'].drop(columns=['variant_type'])

    # Filter somatic_tumour: only keep INS if TYPE is PRECISE
    somatic_tumour = somatic_tumour[
        (somatic_tumour['SVTYPE'] != 'INS') | 
        ((somatic_tumour['SVTYPE'] == 'INS') & (somatic_tumour['TYPE'] == 'PRECISE'))
    ]
    # Create germline_normal dataframe from normal
    germline_normal = normal.dropna(subset=['variant_ID'])
    if 'variant_type' in germline_normal.columns:
        germline_normal = germline_normal.drop(columns=['variant_type'])

    # Identify and create other_normal dataframe
    other_normal = normal[normal['variant_ID'].isna()]
    if 'variant_ID' in other_normal.columns:
        other_normal = other_normal.drop(columns=['variant_ID'])
    if 'variant_type' in other_normal.columns:
        other_normal = other_normal.drop(columns=['variant_type'])

    # Restore original pandas setting
    pd.options.mode.chained_assignment = original_setting

    return somatic_tumour, germline_tumour, germline_normal, other_normal

