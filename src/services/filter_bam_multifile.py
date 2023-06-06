#######################################################################################
# This script filters a BAM file to only include reads that are in a list of read ids #
# Usage: python3 read-reader.py <input.bam> <read list> [suffix]                      #
# If a suffix is provided, the output file will be named <input.filtered.<suffix>.bam #
# Created by andrewprzh                                                               #
# modified to produce multiple filtered bam-files in one run by heidi-holappa         #
#######################################################################################

# TODO: numpy is not used, remove it
import os
import sys
import pysam
import argparse
import numpy
from pathlib import Path


# TODO: Remove this function after verifying that the new function works
def read_reads(inf):
    read_set = set()
    for l in open(inf):
        read_set.add(l.strip().split('\t')[0])
    return read_set

# TODO: Remove this function after verifying that the new function works
def read_called_barcodes(inf):
    read_set = set()
    read_dict = {}
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) >= 7 and v[6] == "*":
            continue
        read_set.add(v[0])
        read_dict[v[0]] = v[6]
    print("Loaded %d read ids" % len(read_dict))
    return read_dict

def create_read_dict(output_filename_dict: dict, original_read_list: str):
    """
    Creates a dictionary of read ids as keys and the transcripts they are mapped to.
    The relationship can be one-to-many, i.e. one read can be mapped to multiple transcripts.

    Args:
        output_filename_dict (dict): a dictionary with transcript ids as keys and output filenames as values
        original_read_list (str): a filename of the original read list

    Returns:
        _type_: _description_
    """
    read_dict = {}
    with open(original_read_list) as file:
        for line in file:
            if line.rstrip("\n").split("\t")[1] in output_filename_dict:
                if line.rstrip("\n").split("\t")[0] not in read_dict:
                    read_dict[line.rstrip("\n").split("\t")[0]] = set()
                read_dict[line.rstrip("\n").split("\t")[0]].add(line.rstrip("\n").split("\t")[1])
    return read_dict


def create_output_filename_dict(bam_file: str, transcript_list: str, suffix: str):
    """
    Creates a dictionary of transcript ids as keys and output filenames as values

    Args:
        bam_file (str): filename of the input bam file
        transcript_list (str): filename of the transcript list
        suffix (str): optional suffix for the output filenames

    Returns:
        dict: transcript ids as keys and output filenames as values
    """
    filenames_dict = {}
    transcripts = []
    with open(transcript_list) as file:
        transcripts = [line.rstrip("\n") for line in file]
    for transcript in transcripts:
       filenames_dict[transcript] = os.path.dirname(os.path.abspath(transcript_list)) + "/" + Path(bam_file).stem + "." + Path(transcript).stem + suffix + ".bam"
    return filenames_dict

def filter_reads(in_file_name: str, out_file_dict: dict, read_dict: dict):
    """
    Filters reads from a BAM file to multiple BAM files based on a dictionary of transcript ids and output filenames

    Args:
        in_file_name (str): filname of the input BAM file
        out_file_dict (dict): a dictionary with transcript ids as keys and output filenames as values
        read_set (dict): a dictionary with read ids as keys and sets of transcript ids as values
    """
    inf = pysam.AlignmentFile(in_file_name, "rb")
    out_files = {}
    for transcript, filename in out_file_dict.items():
        out_files[transcript] = pysam.AlignmentFile(filename, "wb", template=inf)

    count = 0
    passed = 0

    for read in inf:
        if read.reference_id == -1 or read.is_secondary:
            continue

        count += 1
        if count % 10000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        if read.query_name in read_dict:
            for file in read_dict[read.query_name]:
                if file in out_files:
                    out_files[file].write(read)
                else:
                    print("Transcript " + file + " not found in output file dictionary")
            passed += 1

    print("Processed " + str(count) + " reads, written " + str(passed))
    inf.close()
    idx_count = 0
    for transcript, filename in out_files.items():
        idx_count += 1
        sys.stdout.write("Indexing file: " + str(idx_count) + "\r")
        filename.close()
        pysam.index(out_file_dict[transcript])


def init_parser():
    parser = argparse.ArgumentParser(
        description='Filter reads from a BAM file to multiple BAM files based on a list of transcript ids',
        usage='python3 read-reader.py <input.bam> <read list> [suffix]'
    )
    parser.add_argument('-i', '--input', help='input BAM file')
    parser.add_argument('-m', '--model_reads_tsv', help='tsv-file with read ids and transcript ids')
    parser.add_argument('-t', '--transcript_list', help='transcript list. One transcript id per line')
    parser.add_argument('-s', '--suffix', nargs='?', help='optional suffix for the output filenames')
    return parser

if __name__ == "__main__":
    parser = init_parser()

    args = parser.parse_args()

    out_file_dict = create_output_filename_dict(args.input, args.transcript_list, args.suffix)
    read_dict = create_read_dict(out_file_dict, args.model_reads_tsv)
    filter_reads(args.input, out_file_dict, read_dict)