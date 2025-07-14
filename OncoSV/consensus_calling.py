#!/usr/bin/env python3

import pandas as pd
import numpy as np

def remove_empty_columns(df):
    """Remove columns that are entirely empty or all-NA."""
    return df.dropna(axis=1, how='all')
    
def consensus_calling(sniffle, cutesv, svim, chroms, length=20, sd_threshold=0.2):

    # Save original setting
    original_setting = pd.options.mode.chained_assignment

    # Suppress SettingWithCopyWarning
    pd.options.mode.chained_assignment = None
    
    # Remove empty columns from each individual DataFrame
    sniffle = remove_empty_columns(sniffle)
    cutesv = remove_empty_columns(cutesv)
    svim = remove_empty_columns(svim)

    merged_df = pd.concat([sniffle, cutesv, svim], ignore_index=True)
    merged_df = merged_df[(merged_df['CHROM'].isin(chroms)) & (merged_df['CHROM2'].isin(chroms))]

    # Consensus_ID and flag columns
    merged_df['ConsensusSV_ID'] = None
    merged_df['flag'] = 0

    n = 0  # Resetting SV number

    # Process each chromosome
    for chr in chroms:
        chr_df = merged_df[(merged_df['CHROM'] == chr) | (merged_df['CHROM2'] == chr)]

        # Assign values for sorting and applying the new condition
        chr_df[['chrom_sort', 'chrom2_sort', 'pos_sort', 'end_sort']] = \
            chr_df.apply(lambda row: pd.Series({
                'chrom_sort': row['CHROM'] if row['CHROM'] == chr else row['CHROM2'],
                'chrom2_sort': row['CHROM2'] if row['CHROM'] == chr else row['CHROM'],
                'pos_sort': row['POS'] if row['CHROM'] == chr else row['END'],
                'end_sort': row['END'] if row['CHROM'] == chr else row['POS']
            }), axis=1)

        chr_df = chr_df.sort_values(by='pos_sort')

        # First check for large DEL and DUP based on the new condition
        large_del_dup = chr_df[(chr_df['SVTYPE'].isin(['DEL', 'DUP'])) & (chr_df['SVLEN'].abs() > 1000)]
        for i in large_del_dup.index:
            if chr_df.loc[i, 'flag'] == 0:  # Ensure not already processed
                # Find overlapping variants
                overlapping = chr_df[(chr_df['pos_sort'] >= chr_df.loc[i, 'pos_sort'] - length) & 
                                     (chr_df['pos_sort'] <= chr_df.loc[i, 'pos_sort'] + length) &
                                     (chr_df['SVTYPE'] == chr_df.loc[i, 'SVTYPE'])]
                # Calculate the standard deviation of lengths
                if (overlapping['SVLEN'].std() / abs(chr_df.loc[i, 'SVLEN']) < sd_threshold):
                    # Assign the same consensus ID if the condition is met
                    consensus_id = f"consensusSV.type.{n}"
                    for idx in overlapping.index:
                        chr_df.at[idx, 'ConsensusSV_ID'] = consensus_id
                        chr_df.at[idx, 'flag'] = 1
                n += 1

        # Process other variants not yet flagged
        for i in chr_df.index:
            if chr_df.loc[i, 'flag'] == 0:
                n += 1
                chr_df.at[i, 'ConsensusSV_ID'] = f"consensusSV.type.{n}"
                chr_df.at[i, 'flag'] = 1

                # Overlapping check for all variants
                overlapping_rows = chr_df[(chr_df['pos_sort'].isin(range(chr_df.loc[i, 'pos_sort'] - length, chr_df.loc[i, 'pos_sort'] + length))) &
                                          (chr_df['chrom2_sort'] == chr_df.loc[i, 'chrom2_sort']) &
                                          (chr_df['end_sort'].isin(range(chr_df.loc[i, 'end_sort'] - length, chr_df.loc[i, 'end_sort'] + length)))]
                
                chr_df.loc[overlapping_rows.index, 'ConsensusSV_ID'] = f"consensusSV.type.{n}"
                chr_df.loc[overlapping_rows.index, 'flag'] = 1

        # Update the main DataFrame with the results from this chromosome
        merged_df.loc[chr_df.index, 'ConsensusSV_ID'] = chr_df['ConsensusSV_ID']
        merged_df.loc[chr_df.index, 'flag'] = chr_df['flag']

    # Remove unnecessary columns
    merged_df = merged_df.drop(columns=['flag'])
    id_column = merged_df.pop('ConsensusSV_ID')
    merged_df.insert(5, 'ConsensusSV_ID', id_column)

    return merged_df

# Usage:
# result = consensus_calling(df_sniffle, df_cutesv, df_svim, chroms)
