import os

from pathlib import Path
import argparse
import json

from config import SAVE_DIR
from config import DEFAULT_WINDOW_SIZE
from config import CREATE_IMG_N_TRESHOLD


def init_argparser():
    parser = argparse.ArgumentParser(
        description='compAna: a tool for comparing annotations',
        usage='python3 src/compana.py -g <gffcompare_gtf> -r <reference_gtf> -a <reference_fasta> [-o <offset>] [-f] [-s] [-c <class_code>] [-j <json>] [-b <reads_bam>] [-t <reads_tsv>] [-w <window_size>] [-e]'
    )
    parser.add_argument(
        '-g', '--gffcompare_gtf',
        help='GTF-file to be imported into the database',
        required=False,
        metavar='')
    parser.add_argument(
        '-r', '--reference_gtf',
        help='reference GTF-file to be compared against',
        required=False,
        metavar='')
    parser.add_argument(
        '-i', '--isoquant_gtf',
        help="IsoQuant GTF-file to be imported into database",
        required=False,
        metavar=''
    )
    parser.add_argument(
        '-a', '--reference_fasta',
        help='reference FASTA-file to be compared against',
        required=False,
        metavar='')
    parser.add_argument(
        '-o', '--offset',
        help='offset from which canonical splice point is to be searched',
        nargs='+',
        metavar='')
    parser.add_argument(
        '-f', '--force',
        help='force overwrite of existing database',
        action='store_true')
    parser.add_argument(
        '-s', '--stats',
        help='output statistics of class codes',
        action='store_true')
    parser.add_argument(
        '-c', '--class-code',
        nargs='+',
        help='specify gffcompare class code to analyze.')
    parser.add_argument(
        '-j', '--json',
        help='input arguments from json file',
        metavar='')
    parser.add_argument(
        '-b', '--reads_bam',
        help='BAM-file from which reads are to be extracted from',
        metavar='')
    parser.add_argument(
        '-t', '--reads_tsv',
        help='tsv-file for read mapping created by IsoQuant',
        metavar='')
    parser.add_argument(
        '-w', '--window_size',
        help='window from which indels and mismatches are to be searched in interesting locations',
        metavar='',
        default=DEFAULT_WINDOW_SIZE)
    parser.add_argument(
        '-e', '--extended_debug',
        help='enable extended debug output',
        action='store_true')
    parser.add_argument(
        '-m', '--min_reads_for_graph',
        help='threshold for the n of cases for creating images',
        default=CREATE_IMG_N_TRESHOLD)
    parser.add_argument(
        '-n', '--no_canonicals',
        help='canonical splice sites are not considered. Enables more aggressive error prediction and correction.',
        action='store_true')
    parser.add_argument(
        '-v', '--very_conservative',
        help='canonical splice sites are considered, threshold must exceed, there must be a consentration of deletions',
        action='store_true'
    )

    parser_args = parser.parse_args()
    parser_dict = vars(parser_args)
    force_bool = parser_args.force
    no_canonicals_bool = parser_args.no_canonicals
    very_conservative_bool = parser_args.very_conservative

    if parser_args.json:
        with open(parser_args.json, encoding="UTF-8") as json_file:
            json_dict = json.load(json_file)
            parser_dict.update(json_dict)
    elif parser_args.offset:
        parser_dict["offset"] = parser_args.offset[0].split(" ")

    if not parser_args.gffcompare_gtf or not parser_args.reference_gtf:
        parser.print_help()
        exit(1)

    if parser_args.json:
        parser_dict['json_file'] = Path(parser_args.json).stem
        parser_dict['save_file'] = os.path.join(
            SAVE_DIR, Path(parser_args.json).stem + '-matching-cases.pkl')
        parser_dict['intron_save_file'] = os.path.join(
            SAVE_DIR, Path(parser_args.json).stem + '-intron-cases.pkl')
    else:
        parser_dict['save_file'] = os.path.join(SAVE_DIR, Path(
            parser_args.gffcompare_gtf).stem, 'matching-cases.pkl')
        parser_dict['intron_save_file'] = os.path.join(SAVE_DIR, Path(
            parser_args.gffcompare_gtf).stem, 'intron-cases.pkl')

    if not parser_args.offset:
        parser_dict["offset"] = (0, float('inf'))
    else:
        offset_range = []
        for element in parser_args.offset:
            if isinstance(element, str) and len(element) == 3:
                offset_range.append(float('inf'))
            else:
                offset_range.append(int(element))
        parser_dict["offset"] = (
            max(0, offset_range[0]), max(0, offset_range[1])
        )

    if force_bool:
        parser_dict["force"] = True
    if no_canonicals_bool:
        parser_dict["no_canonicals"] = True
    if very_conservative_bool:
        parser_dict["very_conservative"] = True

    return parser_args
