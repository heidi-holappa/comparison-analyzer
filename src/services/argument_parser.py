import argparse
import json


def init_argparser():
    parser = argparse.ArgumentParser(
        description='compAna: a tool for comparing annotations',
        usage='python3 compAna.py -i <input.gtf> [-f] [-s]'
    )
    parser.add_argument(
        '-g', '--gffcompare_gtf',
        help='GTF-file to be imported into the database',
        required=True,
        metavar='')
    parser.add_argument(
        '-r', '--reference_gtf',
        help='reference GTF-file to be compared against',
        required=True,
        metavar='')
    parser.add_argument(
        '-a', '--reference_fasta',
        help='reference FASTA-file to be compared against',
        required=False,
        metavar='')
    parser.add_argument(
        '-o', '--offset',
        help='offset from which canonical splice point is to be searched',
        required=False,
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

    parser = parser.parse_args()
    parser_dict = vars(parser)

    if parser.json:
        with open(parser.json, encoding="UTF-8") as json_file:
            json_dict = json.load(json_file)
            parser_dict.update(json_dict)
            print(parser_dict)
    return parser
