#!/usr/bin/env python3

import gzip

# Function to read VCF lines based on the prefix and optional chromosome filtering
def read_vcf_lines(file_path, line_prefix, chroms=None, extended_chroms=False):
    # Use gzip if the file is compressed; otherwise, use regular open
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'

    with open_func(file_path, mode) as file:
        if chroms:
            if extended_chroms:
                lines = [line.strip() for line in file if line.startswith(line_prefix) and any(line.startswith(line_prefix + chrom + ',') for chrom in chroms)]
            else:
                lines = [line.strip() for line in file if line.startswith(line_prefix) and any(chrom + ',' in line for chrom in chroms)]
        else:
            lines = [line.strip() for line in file if line.startswith(line_prefix)]
    return lines

# Function to combine VCF lines, avoiding duplicates and sorting them
def combine_vcf_lines(file1, file2, line_prefix, chroms=None, extended_chroms=False):
    lines1 = read_vcf_lines(file1, line_prefix, chroms, extended_chroms)
    lines2 = read_vcf_lines(file2, line_prefix, chroms, extended_chroms)

    # Combine lines while avoiding duplicates
    combined_lines = list(set(lines1 + lines2))

    # Sort combined lines if they are contigs
    if chroms:
        primary_chrom_order = {chrom: i for i, chrom in enumerate(chroms)}
        combined_lines.sort(key=lambda x: (primary_chrom_order.get(x.split('ID=')[1].split(',')[0], len(primary_chrom_order)), x))

    return combined_lines
