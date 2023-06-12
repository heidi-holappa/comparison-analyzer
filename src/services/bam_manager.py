import os
from pathlib import Path

from services.filter_bam_multifile import create_output_filename_dict
from services.filter_bam_multifile import create_read_dict
from services.filter_bam_multifile import filter_reads
from services.output_manager import default_output_manager as output_manager
from services.alignment_parser import default_alignment_parser as alignment_parser

from config import TEMPORARY_DIR, CIGAR_RESULTS_LOG


class BamManager:

    def __init__(self, bam_path: str, tsv_path: str, matching_cases_dict: dict):
        self.bam_path = bam_path
        self.tsv_path = tsv_path
        self.matching_cases_dict = matching_cases_dict
        self.transcript_set = set()
        for row in matching_cases_dict:
            self.transcript_set.add(row[0])

    def execute(self):
        output_manager.output_line({
            "line": "PROCESSING BAM-FILE",
            "is_title": True
        })
        output_manager.output_line({
            "line": f"Input BAM-file: {self.bam_path}",
            "is_info": True
        })

        output_filename_dict = create_output_filename_dict(
            self.bam_path,
            self.transcript_set,
            TEMPORARY_DIR
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
        self.remove_temporary_path()
        # TODO: open each file and extract the reads
        # TODO: compare coordinates of reads with coordinates of matching_cases_dict
        # TODO: if there is an indel at the given position, do something. Perhaps calculate percentage of reads with indel?

    def create_temporary_path(self):
        if not os.path.exists(TEMPORARY_DIR):
            os.mkdir(TEMPORARY_DIR)

    def remove_temporary_path(self):
        if os.path.exists(TEMPORARY_DIR):
            for file in os.listdir(TEMPORARY_DIR):
                os.remove(os.path.join(TEMPORARY_DIR, file))
            os.rmdir(TEMPORARY_DIR)

    def iterate_extracted_files(self):
        for key, value in self.matching_cases_dict.items():
            filename = Path(self.bam_path).stem + "." + key[0] + ".bam"
            print(filename, value)
            print(os.listdir(TEMPORARY_DIR))
            print(filename in os.listdir(TEMPORARY_DIR))
            if filename in os.listdir(TEMPORARY_DIR):
                alignment_parser.execute(filename, location=value)

        output_manager.output_line({
            "line": "Insertions and deletions found at given locations",
            "is_info": True
        })

        for key, value in alignment_parser.case_count.items():
            output_manager.output_line({
                "line": f"{key}: {value}",
                "is_info": True
            })
