#!/usr/bin/env python3
import pandas as pd
import re
from collections import defaultdict

def convert_to_list(x):
    if isinstance(x, str):
        try:
            return eval(x)
        except:
            return []  # Handle conversion failure gracefully
    return x

def find_shared_reads(data):
    shared_reads = {}
    for index, row in data.iterrows():
        sv_id = row['ID']
        consensus_id = row['ConsensusSV_ID']
        reads = convert_to_list(row['RNAMES'])
        chrom = row['CHROM']
        chrom2 = row['CHROM2']
        pos = row['POS']
        end = row['END']
        af = row['AF']
        sample = row['Sample']
        for read in reads:
            if read not in shared_reads:
                shared_reads[read] = {
                    'SV_IDs': [],
                    'Consensus_IDs': [],
                    'CHROM_POS': {},  # Stores POS per CHROM
                    'CHROM2_END': {}, # Stores END per CHROM2
                    'CHROM_POS_END': {},  # Stores (pos, end, sv_id) per chromosome
                    'AF': [],  # Store 'AF' as a list
                    'Samples': []  
                }
            shared_reads[read]['SV_IDs'].append(sv_id)
            shared_reads[read]['Consensus_IDs'].append(consensus_id)
            shared_reads[read]['CHROM_POS'].setdefault(chrom, []).append(pos)
            shared_reads[read]['CHROM2_END'].setdefault(chrom2, []).append(end)
            shared_reads[read]['CHROM_POS_END'].setdefault(chrom, []).append((pos, end, sv_id))
            shared_reads[read]['AF'].append(af)
            shared_reads[read]['Samples'].append(sample)
    return shared_reads

def aggregate_chromosomes(chroms):
    return ";".join(sorted(chroms))

def aggregate_positions_by_chromosome(pos_dict):
    aggregated = {chrom: f"{min(positions)}-{max(positions)}" for chrom, positions in pos_dict.items()}
    return ";".join([f"{k}:{v}" for k, v in aggregated.items()])

def extract_sv_type(sv_id):
    # Extract the SV type from the ID, assuming it's always present as a substring
    for sv_type in ['INS', 'DEL', 'DUP', 'INV', 'BND']:
        if sv_type in sv_id:
            return sv_type
    return 'UNK'  # Return 'UNK' for unknown SV types

def aggregate_all_positions(pos_end_dict):
    """ Aggregate all positions into a single string separated by semicolons for each chromosome. """
    pos_aggregated = []
    end_aggregated = []
    sv_count = {}  # Dictionary to keep track of counts for each SV type

    for chrom, positions in pos_end_dict.items():
        for pos, end, sv_id in positions:
            sv_type = extract_sv_type(sv_id)
            if sv_type not in sv_count:
                sv_count[sv_type] = 0
            sv_count[sv_type] += 1
            sv_type_with_count = f"{sv_type}{sv_count[sv_type]}"
            pos_aggregated.append(f"{sv_type_with_count}-{chrom}:{pos}")
            end_aggregated.append(f"{sv_type_with_count}-{chrom}:{end}")

    return ";".join(pos_aggregated), ";".join(end_aggregated)

def format_af_value(af):
    try:
        # Check if 'af' is a tuple and try to unpack and format the first element.
        if isinstance(af, tuple) and af:
            value = float(af[0])  # Convert to float to handle cases where it might not be a proper float.
            return f"{value:.4f}"
        elif isinstance(af, float):  # Direct float handling
            return f"{af:.4f}"
        elif isinstance(af, str) and af.replace('.', '', 1).isdigit():  # Check if 'af' is a numeric string.
            return f"{float(af):.4f}"
    except (TypeError, ValueError) as e:
        # Handle cases where conversion is not possible or 'af' is None.
        return "0.0000"  # Provide a default format if 'af' is incorrect.
    return "0.0000" 

def process_shared_reads(dataframe):
    dataframe['RNAMES'] = dataframe['RNAMES'].apply(convert_to_list)
    shared_read_mapping = find_shared_reads(dataframe)
    output_data = []
    for read, details in shared_read_mapping.items():
        if len(details['SV_IDs']) > 1:  # Only consider shared reads
            pos_bkpt, end_bkpt = aggregate_all_positions(details['CHROM_POS_END'])
            af_values = ';'.join(format_af_value(af) for af in details['AF'])

            output_data.append({
                'Read': read,
                'ID': ','.join(set(details['SV_IDs'])),
                'ConsensusSV_ID': ','.join(set(details['Consensus_IDs'])),
                'CHROM': aggregate_chromosomes(details['CHROM_POS'].keys()),
                'CHROM2': aggregate_chromosomes(details['CHROM2_END'].keys()),
                'POS': aggregate_positions_by_chromosome(details['CHROM_POS']),
                'END': aggregate_positions_by_chromosome(details['CHROM2_END']),
                'POS_BKPT': pos_bkpt,
                'END_BKPT': end_bkpt,
                'AF': af_values,
                'Sample': ",".join(set(details['Samples']))
            })
    return pd.DataFrame(output_data)

def summarize_sv_types(df):
    def get_csv_type(sv_ids):
        types_count = {}
        for sv_id in sv_ids.split(','):
            sv_type = extract_sv_type(sv_id)
            if sv_type in types_count:
                types_count[sv_type] += 1
            else:
                types_count[sv_type] = 1

        summary = []
        for sv_type, count in types_count.items():
            if count > 1:
                summary.append(f"({count}){sv_type}")
            else:
                summary.append(sv_type)

        return '+'.join(summary)

    df['CSV_Type'] = df['ID'].apply(get_csv_type)
    return df

def process_sv_data_with_sv_count(data):
    # Group by 'ID' and count the number of reads
    read_counts = data.groupby('ID')['Read'].count().reset_index(name='Read_Count')
    # Count the number of SVs in each combination
    read_counts['SV_Count'] = read_counts['ID'].apply(lambda x: len(x.split(',')))

    # Assuming details to be merged based on unique 'ID', aggregating 'POS' and 'END' per 'CHROM' and 'CHROM2'
    details = data[['ID', 'CHROM', 'CHROM2', 'POS', 'END', 'POS_BKPT', 'END_BKPT', 'AF']].drop_duplicates('ID')

    # Merge these details back into the read_counts DataFrame
    read_counts = read_counts.merge(details, on='ID', how='left')

    # Ensure that 'Sample' is present in the original 'data' DataFrame before using it
    if 'Sample' in data.columns:
        read_counts['Sample'] = data['Sample']

    columns_order = ['ID', 'CHROM', 'CHROM2', 'POS', 'END', 'POS_BKPT', 'END_BKPT', 'Read_Count', 'SV_Count', 'AF']
    if 'Sample' in read_counts.columns:
        columns_order.append('Sample')

    read_counts = read_counts[columns_order]

    # Add the CSV_Type column
    read_counts = summarize_sv_types(read_counts)

    return read_counts

def extract_chr_and_coordinate(bp):
    match = re.search(r'(chr[\dXY]+):(\d+)', bp)
    return (match.group(1), int(match.group(2))) if match else (None, None)

def process_breakpoints(df):
    def process_row(pos_bkpt, end_bkpt):
        pos_list = [f"{bp.split('-')[0]}-pos-{bp.split('-')[1]}" for bp in pos_bkpt.split(';')]
        end_list = [f"{bp.split('-')[0]}-end-{bp.split('-')[1]}" for bp in end_bkpt.split(';')]

        pos_dict = defaultdict(list)
        end_dict = defaultdict(list)

        for bp in pos_list:
            chr, coord = extract_chr_and_coordinate(bp)
            if chr:
                pos_dict[chr].append(bp)

        for bp in end_list:
            chr, coord = extract_chr_and_coordinate(bp)
            if chr:
                end_dict[chr].append(bp)

        final_order = []

        for chr in sorted(pos_dict.keys() | end_dict.keys()):
            combined_list = pos_dict[chr] + end_dict[chr]
            sorted_list = sorted(combined_list, key=lambda x: extract_chr_and_coordinate(x)[1])
            final_order.append("_".join(sorted_list))

        return ";".join(final_order)

    df['final_combination'] = df.apply(lambda row: process_row(row['POS_BKPT'], row['END_BKPT']), axis=1)
    return df

def check_overlapping_sv(final_combination):
    chr_combinations = final_combination.split(';')
    overlap_results = []

    for combination in chr_combinations:
        breakpoints = combination.split('_')
        open_sv = []
        overlaps = []

        for bp in breakpoints:
            sv = bp.split('-')[0]
            if 'pos' in bp:
                open_sv.append(sv)
            elif 'end' in bp:
                if sv in open_sv:
                    open_sv.remove(sv)
            if len(open_sv) > 1:
                for sv in open_sv:
                    if sv not in overlaps:
                        overlaps.append(sv)

        if overlaps:
            overlap_results.append("-".join(overlaps))
        else:
            overlap_results.append("no")

    return ";".join(overlap_results)

def add_overlapping_column(df):
    df['any_overlapping_sv'] = df['final_combination'].apply(check_overlapping_sv)
    return df


