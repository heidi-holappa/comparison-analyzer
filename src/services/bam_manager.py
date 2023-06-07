import os
from pathlib import Path

from services.filter_bam_multifile import create_output_filename_dict
from services.filter_bam_multifile import create_read_dict
from services.filter_bam_multifile import filter_reads
from services.output_manager import default_output_manager as output_manager
from services.cigar_parser import CigarParser


class BamManager:

    def __init__(self, bam_path: str, tsv_path: str, matching_cases_dict: dict):
        self.bam_path = bam_path
        self.tsv_path = tsv_path
        self.matching_cases_dict = matching_cases_dict
        self.transcript_set = set()
        for row in matching_cases_dict:
            self.transcript_set.add(row[0])
        self.temporary_path = self.create_temporary_path()

    def execute(self):
        output_manager.output_line("PROCESSING BAM-FILE", is_title=True)
        output_manager.output_line(
            "Fetching reads from BAM-file. This might take some time.", is_info=True)
        output_filename_dict = create_output_filename_dict(
            self.bam_path,
            self.transcript_set,
            self.temporary_path
        )
        read_dict = create_read_dict(
            output_filename_dict,
            self.tsv_path
        )
        filter_reads(
            self.bam_path,
            output_filename_dict,
            read_dict
        )
        self.iterate_extracted_files()
        os.rmdir(self.temporary_path)
        # TODO: open each file and extract the reads
        # TODO: compare coordinates of reads with coordinates of matching_cases_dict
        # TODO: if there is an indel at the given position, do something. Perhaps calculate percentage of reads with indel?
        # TODO: delete temporary files

    def create_temporary_path(self):
        temporary_dir = "temporary_files"
        current_dir = os.getcwd()
        temporary_path = os.path.join(current_dir, temporary_dir)
        if not os.path.exists(temporary_path):
            os.makedirs(temporary_path)
        return temporary_path

    def iterate_extracted_files(self):
        cigar_parser = CigarParser()
        cigar_results = {}
        for key, value in self.matching_cases_dict.items():
            filename = Path(self.bam_path).stem + key[0] + ".bam"
            if filename in os.listdir(self.temporary_path):
                samfile = cigar_parser.initialize_file(filename)
                cigar_results[filename] = cigar_parser.extract_cigar_symbol(
                    samfile, value)

        with open("cigar_results.txt", "w") as file:
            for key, value in cigar_results.items():
                file.write(f"{key}: {value}\n")
