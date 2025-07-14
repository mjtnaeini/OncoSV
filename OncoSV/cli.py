#!/usr/bin/env python3

import argparse
from .main_consensus import run_consensus
from .main_somatic import run_pair
from .main_complexSV import run_complexSV

def main():
    parser = argparse.ArgumentParser(description="ComplexSVnet Package CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subparser for the 'consensus' command
    parser_consensus = subparsers.add_parser(
        'consensus',
        help='Run consensus structural variant calling for individual samples'
    )
    parser_consensus.add_argument('-s', '--sniffles', type=str, required=True, help='Sniffles VCF file')
    parser_consensus.add_argument('-c', '--cutesv', type=str, required=True, help='CuteSV VCF file')
    parser_consensus.add_argument('-v', '--svim', type=str, required=True, help='SVIM VCF file')
    parser_consensus.add_argument('-o', '--out-file', type=str, required=True, help='Output VCF file')

    # Optional arguments for 'consensus'
    parser_consensus.add_argument('-x', '--chrom', type=str, help='Which chromosomes to query', default='all')
    parser_consensus.add_argument('--sample-id', type=str, help='Sample ID', default='Sample')
    parser_consensus.add_argument('-q', '--quality-threshold', type=int, help='Minimum quality of SVs', default=10)
    parser_consensus.add_argument('-m', '--minimum-sv-size', type=int, help='Minimum SV size', default=50)
    parser_consensus.add_argument('-M', '--maximum-sv-size', type=int, help='Maximum SV size', default=1000000)
    parser_consensus.add_argument('--compress', action='store_true', help='Compress the VCF file')
    parser_consensus.add_argument('--apply-af-filtering', type=str, choices=["true", "false"], help='AF filtering')

    # Subparser for the 'pair' command
    parser_pair = subparsers.add_parser('pair', help='Run somatic and germline variant calling for paired samples')

    # Required arguments for 'pair'
    parser_pair.add_argument('-t', '--tumour-consensus', type=str, required=True, help='Tumour VCF file')
    parser_pair.add_argument('-o', '--out-dir', type=str, required=True, help='Output directory')

    # New argument for normal mode selection
    parser_pair.add_argument(
        '--normal-mode', type=str, choices=['single', 'multi'], required=True,
        help="Specify if using a single ('single') or multiple ('multi') normal sample VCFs"
    )

    # Arguments for normal sample VCF files
    parser_pair.add_argument('-n', '--normal-sample', type=str, help='Single normal sample VCF (if normal-mode=single)')
    parser_pair.add_argument('--normal-sample1', type=str, help='First normal sample VCF (if normal-mode=multi)')
    parser_pair.add_argument('--normal-sample2', type=str, help='Second normal sample VCF (if normal-mode=multi)')
    parser_pair.add_argument('--normal-sample3', type=str, help='Third normal sample VCF (if normal-mode=multi)')
    
    # Option to save merged normal samples
    parser_pair.add_argument('--save-merged-normal', type=str, choices=["true", "false"], default="false", help="Save merged normal samples as CSV (default: false)")

    # Optional arguments for 'pair'
    parser_pair.add_argument('-x', '--chrom', type=str, help='Chromosomes to query', default='all')
    parser_pair.add_argument('--vcf-format', type=str, help='VCF format', default='consensus')
    parser_pair.add_argument('--tumour-id', type=str, help='Tumour sample ID', default='Sample')
    parser_pair.add_argument('--normal-id', type=str, help='Normal sample ID', default='Sample')
    parser_pair.add_argument('-q', '--quality-threshold', type=int, help='Minimum quality of SVs', default=10)
    parser_pair.add_argument('-sv', '--minimum-sv-size', type=int, help='Minimum SV size', default=50)
    parser_pair.add_argument('-M', '--maximum-sv-size', type=int, help='Maximum SV size', default=1000000)
    parser_pair.add_argument('--only-somatic', action='store_true', help='Only generate VCF for somatic variants')
    parser_pair.add_argument('--compress', action='store_true', help='Compress VCF file')
    parser_pair.add_argument('--patient-id', type=str, help='Patient ID label')
    parser_pair.add_argument('--svcaller', type=str, choices=["consensus", "sniffles", "cutesv", "svim"], help="SV caller", default="consensus")

    # Subparser for the 'complexSV' command
    parser_complexSV = subparsers.add_parser('complexSV', help='Run complex SV analysis')

    # Required arguments for 'complexSV'
    parser_complexSV.add_argument("--vcf", type=str, required=True, help="Path to the VCF file")
    parser_complexSV.add_argument("--output_dir", type=str, required=True, help="Output directory")

    # Optional arguments for 'complexSV'
    parser_complexSV.add_argument('-x', '--chrom', type=str, help='Chromosomes to query', default='all')
    parser_complexSV.add_argument("--sample_id", type=str, help="Sample ID")
    parser_complexSV.add_argument("--qual", type=int, default=10, help="Quality threshold")
    parser_complexSV.add_argument('-sv', '--minimum-sv-size', type=int, help='Minimum SV size', default=50)
    parser_complexSV.add_argument('-M', '--maximum-sv-size', type=int, help='Maximum SV size', default=1000000)
    parser_complexSV.add_argument("--vcf_format", default='consensus', help="VCF file format")
    parser_complexSV.add_argument("--label_prefix", type=str, default='', help="Label prefix for output filenames")
    parser_complexSV.add_argument('--svcaller', type=str, choices=["consensus", "sniffles", "cutesv", "svim"], help="SV caller", default="consensus")

    args = parser.parse_args()

    if args.command == "consensus":
        run_consensus(args)
    elif args.command == "pair":
        run_pair(args)
    elif args.command == "complexSV":
        run_complexSV(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

