#!/usr/bin/env python3

import pandas as pd
import re
import pysam
import gzip

def format_genotype(genotype):
    if genotype is None:
        return "."
    elif isinstance(genotype, tuple):
        formatted_genotype = [str(g) if g is not None else '.' for g in genotype]
        return "/".join(formatted_genotype)
    return genotype

def clean_tuple_field(field_value):
    if isinstance(field_value, tuple):
        return ",".join(str(item) for item in field_value)
    return field_value

def create_info_field(row, include_variant_ID=False, svcaller="consensus"):
    # Define the keys to include in the info field
    info_keys = ["SVTYPE", "SVLEN", "END", "RNAMES", "AF", "NUM_CALLERS"]

    # Insert 'ConsensusSV_ID' at the beginning if svcaller is 'consensus'
    if svcaller == "consensus":
        info_keys.insert(0, "ConsensusSV_ID")

    info_parts = []

    # Handle precision flags based on 'TYPE'
    if 'TYPE' in row and pd.notna(row['TYPE']):
        if row['TYPE'].upper() == "PRECISE":
            info_parts.append("PRECISE")
        elif row['TYPE'].upper() == "IMPRECISE":
            info_parts.append("IMPRECISE")

    # Handle inclusion of 'Variant_ID'
    if include_variant_ID and 'variant_ID' in row and not pd.isna(row['variant_ID']):
        if svcaller == "consensus":
            idx = info_keys.index("ConsensusSV_ID") + 1
            info_keys.insert(idx, "variant_ID")
        else:
            info_keys.insert(0, "variant_ID")

    for key in info_keys:
        if key == 'variant_ID':
            info_parts.append(f"Variant_ID={row['variant_ID']}")
            continue

        if key in row and not pd.isna(row[key]):
            if key == 'SVLEN' and row.get("SVTYPE") == "BND":
                continue  # Skip SVLEN for BND type SVs
            value = row[key]
            if key == 'NUM_CALLERS':
                value = row.get("NUM_CALLERS", "1")
            elif key == 'SVLEN':
                try:
                    value = int(value)
                except ValueError:
                    value = "."
            elif key == 'AF':
                try:
                    value = float(value)
                    value = f"{value:.2f}"  # Format to two decimal places
                except ValueError:
                    value = "."

            value_str = clean_tuple_field(value) if isinstance(value, tuple) else str(value)
            info_parts.append(f"{key}={value_str}")
        else:
            if key == 'NUM_CALLERS':
                info_parts.append(f"{key}=1")
            elif key == 'SVLEN' and row.get("SVTYPE") != "BND":
                info_parts.append(f"{key}=.")
            else:
                continue

    return ";".join(info_parts)


def create_format_and_sample_fields(row):
    format_keys = ["GT", "GQ", "DR", "DV"]
    format_values = [
        format_genotype(row.get("Genotype", ".")),
        str(int(row.get("GenotypeQuality", "."))) if pd.notna(row.get("GenotypeQuality")) and row.get("GenotypeQuality") != "." else ".",
        str(int(row.get("ReferenceReads", "."))) if pd.notna(row.get("ReferenceReads")) and row.get("ReferenceReads") != "." else ".",
        str(int(row.get("VariantReads", "."))) if pd.notna(row.get("VariantReads")) and row.get("VariantReads") != "." else "."

    ]
    return ":".join(format_keys), ":".join(str(x) for x in format_values)

def generate_vcf_from_dataframe(dataframe, combined_contigs, combined_filters, output_filename, is_compressed=False, sample_id=None):
    # Extract the sample ID from the first row
    if not sample_id:
        sample_id = dataframe.iloc[0]['Sample'] if 'Sample' in dataframe.columns and len(dataframe) > 0 else 'DefaultSample'

    vcf_header = [
        "##fileformat=VCFv4.2",
        '##source=ComplexSVnet-v1.0'
    ]

    # Add contigs to the header
    vcf_header.extend(combined_contigs)
    # Add filters to the header 
    vcf_header.extend(combined_filters)

    # Add specific ALT fields
    vcf_header.extend([
        "##ALT=<ID=INS,Description=\"Insertion\">",
        "##ALT=<ID=DEL,Description=\"Deletion\">",
        "##ALT=<ID=DUP,Description=\"Duplication\">",
        "##ALT=<ID=INV,Description=\"Inversion\">",
        "##ALT=<ID=BND,Description=\"Breakend; Translocation\">"
    ])
    
    vcf_header.extend([
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
        "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Reference Reads\">",
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Variant Reads\">",
    ])

    # Specific INFO fields for consensus
    vcf_header.extend([
        "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise Structural variant\">",
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">",
        "##INFO=<ID=ConsensusSV_ID,Number=1,Type=String,Description=\"Identifier for consensus structural variant\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of structural variant\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
        "##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Names of reads supporting SVs\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
        "##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,Description=\"Number of SV callers reporting this variant\">"      
    ])

    vcf_header.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}")

    mode = 'wb' if is_compressed else 'w'
    open_func = pysam.BGZFile if is_compressed else open

    with open_func(output_filename, mode) as file:
        for line in vcf_header:
            file.write(f"{line}\n".encode('utf-8') if is_compressed else f"{line}\n")
        for _, row in dataframe.iterrows():
            info_field = create_info_field(row)
            format_field, sample_field = create_format_and_sample_fields(row)
            row["QUAL"] = round(float(row["QUAL"])) if pd.notna(row["QUAL"]) else "."
            row["ALT"] = clean_tuple_field(row["ALT"])
            vcf_row = [
                row["CHROM"],
                row["POS"],
                row["ID"],
                row["REF"],
                row["ALT"],
                row["QUAL"],
                row["FILTER"],
                info_field,
                format_field,
                sample_field
            ]
            file.write("\t".join(map(str, vcf_row)).encode('utf-8') if is_compressed else "\t".join(map(str, vcf_row)) + "\n")

    if is_compressed:
        # Create tabix index
        pysam.tabix_index(output_filename, preset="vcf")

def retrieve_vcf_header(file_path):
    contigs, filters = [], []
    with (gzip.open if file_path.endswith('.gz') else open)(file_path, 'rt', encoding='utf-8') as file:
        for line in file:
            if line.startswith("##contig="):
                contigs.append(line.strip())
            elif line.startswith("##FILTER="):
                filters.append(line.strip())
    return contigs, filters

def generate_vcf_variants(dataframe, header_file_path, output_filename, is_compressed=False, include_variant_ID=True, sample_id=None, svcaller="consensus"):
    if not sample_id:
        sample_id = dataframe.iloc[0]['Sample'] if 'Sample' in dataframe.columns else 'DefaultSample'

    contigs, filters = retrieve_vcf_header(header_file_path)

    vcf_header = [
        "##fileformat=VCFv4.2",
        "##source=ComplexSVnet-v1.0"
    ]
    vcf_header.extend(contigs)
    vcf_header.extend(filters)

    # Add ALT fields based on svcaller
    if svcaller == "svim":
        alt_fields = [
            "##ALT=<ID=DEL,Description=\"Deletion\">",
            "##ALT=<ID=INV,Description=\"Inversion\">",
            "##ALT=<ID=DUP,Description=\"Duplication\">",
            "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">",
            "##ALT=<ID=DUP:INT,Description=\"Interspersed Duplication\">",
            "##ALT=<ID=INS,Description=\"Insertion\">",
            "##ALT=<ID=BND,Description=\"Breakend\">"
        ]
    else:
        alt_fields = [
            "##ALT=<ID=INS,Description=\"Insertion\">",
            "##ALT=<ID=DEL,Description=\"Deletion\">",
            "##ALT=<ID=DUP,Description=\"Duplication\">",
            "##ALT=<ID=INV,Description=\"Inversion\">",
            "##ALT=<ID=BND,Description=\"Breakend; Translocation\">"
        ]
    vcf_header.extend(alt_fields)

    # Dynamically include PRECISE and IMPRECISE in the header if used
    if dataframe["TYPE"].str.upper().isin(["PRECISE", "IMPRECISE"]).any():
        vcf_header.extend([
            "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise Structural variant\">",
            "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">"
        ])

    # Add INFO fields
    info_fields = [
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of structural variant\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
        "##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Names of reads supporting SVs\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
        "##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,Description=\"Number of SV callers reporting this variant\">"
    ]
    
    # Add ConsensusSV_ID if svcaller is consensus
    if svcaller == "consensus":
        info_fields.insert(0, "##INFO=<ID=ConsensusSV_ID,Number=1,Type=String,Description=\"Identifier for consensus structural variant\">")

    # Add Variant_ID after PRECISE/IMPRECISE or after ConsensusSV_ID for consensus
    if "variant_ID" in dataframe.columns:
        variant_id_field = "##INFO=<ID=Variant_ID,Number=1,Type=String,Description=\"Identifier including both somatic and germline variants\">"
        if svcaller == "consensus":
            # Add Variant_ID after ConsensusSV_ID
            idx = next((i for i, field in enumerate(info_fields) if "ConsensusSV_ID" in field), -1)
            info_fields.insert(idx + 1, variant_id_field)
        else:
            # Add Variant_ID after PRECISE and IMPRECISE
            precise_idx = len(vcf_header) + 1
            vcf_header.insert(precise_idx, variant_id_field)

    vcf_header.extend(info_fields)

    # Add FORMAT and mandatory column header
    vcf_header.extend([
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
        "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Reference Reads\">",
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Variant Reads\">",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}"
    ])

    mode = 'wb' if is_compressed else 'w'
    open_func = pysam.BGZFile if is_compressed else open

    with open_func(output_filename, mode) as file:
        for line in vcf_header:
            file.write(f"{line}\n".encode('utf-8') if is_compressed else f"{line}\n")
        for _, row in dataframe.iterrows():
            info_field = create_info_field(row, include_variant_ID=include_variant_ID, svcaller=svcaller)
            format_field, sample_field = create_format_and_sample_fields(row)
            vcf_row = [
                row["CHROM"],
                row["POS"],
                row["ID"],
                row["REF"],
                row["ALT"],
                round(float(row["QUAL"])) if pd.notna(row["QUAL"]) else ".",
                row["FILTER"],
                info_field,
                format_field,
                sample_field
            ]
            file.write("\t".join(map(str, vcf_row)).encode('utf-8') if is_compressed else "\t".join(map(str, vcf_row)) + "\n")

    if is_compressed:
        pysam.tabix_index(output_filename, preset="vcf")
