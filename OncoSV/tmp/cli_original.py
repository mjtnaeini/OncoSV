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

    # Required arguments for 'consensus'
    parser_consensus.add_argument(
        '-s',
        '--sniffles',
        type=str,
        required=True,
        help='Sniffles VCF file'
    )
    parser_consensus.add_argument(
        '-c',
        '--cutesv',
        type=str,
        required=True,
        help='CuteSV VCF file'
    )
    parser_consensus.add_argument(
        '-v',
        '--svim',
        type=str,
        required=True,
        help='SVIM VCF file'
    )
    parser_consensus.add_argument(
        '-o',
        '--out-file',
        type=str,
        required=True,
        help='Output VCF file'
    )

    # Optional arguments for 'consensus'
    parser_consensus.add_argument(
        '-x',
        '--chrom',
        type=str,
        help='Which chromosomes to query. Comma-separated list or [all]',
        default='all'
    )
    parser_consensus.add_argument(
        '--sample-id',
        type=str,
        help='Sample ID',
        default='Sample'
    )
    parser_consensus.add_argument(
        '-q',
        '--quality-threshold',
        type=int,
        help='Minimum quality of structural variants',
        default=10
    )
    parser_consensus.add_argument(
        '-m',
        '--minimum-sv-size',
        type=int,
        help='Minimum structural variants size',
        default=50
    )
    parser_consensus.add_argument(
        '-M',
        '--maximum-sv-size',
        type=int,
        help='Maximum structural variants size',
        default=1000000
    )

    parser_consensus.add_argument(
        '--compress',
        action='store_true',
        help='Compress the VCF file using bgzip and create a .tbi index'
    )

    parser_consensus.add_argument(
    '--apply-af-filtering',
    type=str,
    choices=["true", "false"],
    help='Enable or disable AF filtering (default: true)'
    )

    # Subparser for the 'pair' command
    parser_pair = subparsers.add_parser(
        'pair',
        help='Run somatic and germline variant calling for paired cancer samples'
    )

    # Required arguments for 'pair'
    parser_pair.add_argument(
        '-t',
        '--tumour-consensus',
        type=str,
        required=True,
        help='Path to the tumour sample consensus structural variant VCF file.'
    )
    parser_pair.add_argument(
        '-n',
        '--normal-consensus',
        type=str,
        required=True,
        help='Path to the normal sample consensus structural variant VCF file.'
    )
    parser_pair.add_argument(
        '-o',
        '--out-dir',
        type=str,
        required=True,
        help='Output directory for VCF file'
    )

    # Optional arguments for 'pair'
    parser_pair.add_argument(
        '-x',
        '--chrom',
        type=str,
        required=False,
        help='Which chromosomes to query (comma-separated list or "all")',
        default='all'
    )
    parser_pair.add_argument(
        '--vcf-format',
        type=str,
        required=False,
        help='Format of the VCF files (default is "consensus")',
        default='consensus'
    )
    parser_pair.add_argument(
        '--tumour-id',
        type=str,
        required=False,
        help='Tumour sample ID',
        default='Sample'
    )
    parser_pair.add_argument(
        '--normal-id',
        type=str,
        required=False,
        help='Normal sample ID',
        default='Sample'
    )
    parser_pair.add_argument(
        '-q',
        '--quality-threshold',
        type=int,
        required=False,
        help='Minimum quality of structural variants',
        default=10
    )
    parser_pair.add_argument(
        '-sv',
        '--minimum-sv-size',
        type=int,
        required=False,
        help='Minimum structural variants size',
        default=50
    )
    parser_pair.add_argument(
        '-M',
        '--maximum-sv-size',
        type=int,
        help='Maximum structural variants size',
        default=1000000
    )
    parser_pair.add_argument(
        '--only-somatic',
        action='store_true',
        required=False,
        help='Generate VCF file only for somatic tumour variants'
    )
    parser_pair.add_argument(
        '--compress',
        action='store_true',
        required=False,
        help='Compress the VCF file using bgzip and create a .tbi index'
    )
    parser_pair.add_argument(
        '--patient-id',
        type=str,
        required=False,
        default=None,
        help='Optional patient ID to label in the VCF files'
    )
    parser_pair.add_argument(
        '--svcaller',
        type=str,
        choices=["consensus", "sniffles", "cutesv", "svim"],
        required=False,
        help="Specify the structural variant caller (default is 'consensus')",
        default="consensus"
    )

    # Subparser for the 'complexSV' command
    parser_complexSV = subparsers.add_parser(
        'complexSV',
        help='Run complex structural variant analysis'
    )

    # Required arguments for 'complexSV'
    parser_complexSV.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="Path to the VCF file"
    )
    parser_complexSV.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory for the CSV files"
    )

    # Optional arguments for 'complexSV'
    parser_complexSV.add_argument(
        '-x',
        '--chrom',
        type=str,
        required=False,
        help='Which chromosomes to query. Comma-separated list or [all]',
        default='all'
    )
    parser_complexSV.add_argument(
        "--sample_id",
        type=str,
        required=False,
        help="Sample ID"
    )
    parser_complexSV.add_argument(
        "--qual",
        type=int,
        default=10,
        help="Quality threshold"
    )
    parser_complexSV.add_argument(
        '-sv',
        '--minimum-sv-size',
        type=int,
        required=False,
        help='Minimum structural variants size',
        default=50
    )
    parser_complexSV.add_argument(
        '-M',
        '--maximum-sv-size',
        type=int,
        help='Maximum structural variants size',
        default=1000000
    )

    parser_complexSV.add_argument(
        "--vcf_format",
        default='consensus',
        help="Format of the VCF file"
    )
    parser_complexSV.add_argument(
        "--label_prefix", 
        type=str, 
        default='', 
        help="Optional label prefix for output filenames"
    )
    parser_complexSV.add_argument(
        '--svcaller',
        type=str,
        choices=["consensus", "sniffles", "cutesv", "svim"],
        required=False,
        help="Specify the structural variant caller (default is 'consensus')",
        default="consensus"
    )

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

