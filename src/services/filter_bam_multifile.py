#######################################################################################
# This script filters a BAM file to only include reads that are in a list of read ids #
# Usage: python3 read-reader.py <input.bam> <read list> [suffix]                      #
# If a suffix is provided, the output file will be named <input.filtered.<suffix>.bam #
# Created by andrewprzh                                                               #
# modified to produce multiple filtered bam-files in one run by heidi-holappa         #
#######################################################################################

import os
import argparse
from pathlib import Path

import pysam
from services.output_manager import default_output_manager as output_manager


def create_dict_of_reads(transcript_set: set, original_read_list: str):
    """
    Creates a dictionary of read ids as keys and the transcripts they are mapped to.
    The relationship can be one-to-many, i.e. one read can be mapped to multiple transcripts.

    Args:
        output_filename_dict (dict): a dict with transcript ids as keys, output filenames as values
        original_read_list (str): a filename of the original read list

    Returns:
        dict: outputs a dictionary of read ids as keys and the transcripts they are mapped to
    """
    new_read_dict = {}
    with open(original_read_list, encoding="UTF-8") as file:
        for line in file:
            print(line.rstrip("\n").split("\t")[1])
            if line.rstrip("\n").split("\t")[1] in transcript_set:
                if line.rstrip("\n").split("\t")[0] not in new_read_dict:
                    new_read_dict[line.rstrip("\n").split("\t")[0]] = set()
                new_read_dict[line.rstrip("\n").split("\t")[0]].add(
                    line.rstrip("\n").split("\t")[1])
    return new_read_dict


def create_read_dict(output_filename_dict: dict, original_read_list: str):
    """
    Creates a dictionary of read ids as keys and the transcripts they are mapped to.
    The relationship can be one-to-many, i.e. one read can be mapped to multiple transcripts.

    Args:
        output_filename_dict (dict): a dict with transcript ids as keys, output filenames as values
        original_read_list (str): a filename of the original read list

    Returns:
        dict: outputs a dictionary of read ids as keys and the transcripts they are mapped to
    """
    new_read_dict = {}
    with open(original_read_list, encoding="UTF-8") as file:
        for line in file:
            if line.rstrip("\n").split("\t")[1] in output_filename_dict:
                if line.rstrip("\n").split("\t")[0] not in new_read_dict:
                    new_read_dict[line.rstrip("\n").split("\t")[0]] = set()
                new_read_dict[line.rstrip("\n").split("\t")[0]].add(
                    line.rstrip("\n").split("\t")[1])
    return new_read_dict


def create_output_filename_dict_cli(bam_file: str, transcript_list: str, suffix: str):
    """
    Creates a dictionary of transcript ids as keys and output filenames as values. 
    This method is used when the script is run from the command line.

    Args:
        bam_file (str): filename of the input bam file
        transcript_list (str): filename of the transcript list
        suffix (str): optional suffix for the output filenames

    Returns:
        dict: transcript ids as keys and output filenames as values
    """
    filenames_dict = {}
    transcripts = []
    with open(transcript_list, encoding="UTF-8") as file:
        transcripts = [line.rstrip("\n") for line in file]
    for transcript in transcripts:
        filenames_dict[transcript] = os.path.dirname(os.path.abspath(
            transcript_list)) + "/" + Path(bam_file).stem + "." + \
            transcript + suffix + ".bam"
    return filenames_dict


def create_output_filename_dict(bam_file: str, transcript_set: set, temporary_path: str):
    """
        Creates a dictionary of transcript ids as keys and output filenames as values. 
        This method is used when the script is run as a part of a script.
    """
    filename_dict = {}
    for transcript in transcript_set:
        filename = Path(bam_file).stem + "." + Path(transcript).stem + ".bam"
        filename_dict[transcript] = os.path.join(temporary_path, filename)
    return filename_dict


def filter_reads(bam_file: str, dict_of_output_filenames: dict, dict_of_reads: dict):
    # pylint: disable=no-member
    """
    Filters reads from a BAM file to multiple BAM files 
    based on a dictionary of transcript ids and output filenames

    Args:
        in_file_name (str): filname of the input BAM file
        out_file_dict (dict): a dict with transcript ids as keys and output filenames as values
        read_set (dict): a dictionary with read ids as keys and sets of transcript ids as values
    """
    inf = pysam.AlignmentFile(bam_file, "rb")
    out_files = {}
    for transcript, filename in dict_of_output_filenames.items():
        out_files[transcript] = pysam.AlignmentFile(
            filename, "wb", template=inf)

    count = 0
    passed = 0

    for read in inf:
        if read.reference_id == -1 or read.is_secondary:
            continue

        count += 1
        if count % 10000 == 0:
            output_manager.output_line({
                "line": "Processed " + str(count) + " reads, written " + str(passed),
                "end_line": "\r",
                "is_info": True
            })

        if read.query_name in dict_of_reads:
            for file in dict_of_reads[read.query_name]:
                if file in out_files:
                    out_files[file].write(read)
                else:
                    output_manager.output_line({
                        "line": "Transcript " + file + " not found in output file dictionary",
                        "is_error": True
                    })
            passed += 1

    output_manager.output_line({
        "line": "Processed " + str(count) + " reads, written " + str(passed),
        "is_info": True
    })
    inf.close()
    idx_count = 0
    for transcript, filename in out_files.items():
        idx_count += 1
        output_manager.output_line({
            "line": "Indexing file: " + str(idx_count),
            "end_line": "\r",
            "is_info": True
        })
        filename.close()
        pysam.index(dict_of_output_filenames[transcript])


def init_parser():
    cli_parser = argparse.ArgumentParser(
        description='Filter reads from a BAM file to multiple \
            BAM files based on a list of transcript ids',
        usage='python3 read-reader.py <input.bam> <read list> [suffix]'
    )
    cli_parser.add_argument('-i', '--input', help='input BAM file')
    cli_parser.add_argument('-m', '--model_reads_tsv',
                            help='tsv-file with read ids and transcript ids')
    cli_parser.add_argument('-t', '--transcript_list',
                            help='transcript list. One transcript id per line')
    cli_parser.add_argument('-s', '--suffix', nargs='?',
                            help='optional suffix for the output filenames')
    return cli_parser


if __name__ == "__main__":
    parser = init_parser()

    args = parser.parse_args()

    out_file_dict = create_output_filename_dict_cli(
        args.input, args.transcript_list, args.suffix)
    read_dict = create_read_dict(out_file_dict, args.model_reads_tsv)
    filter_reads(args.input, out_file_dict, read_dict)
