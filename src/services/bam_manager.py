import os
from pathlib import Path

from services.filter_bam_multifile import create_output_filename_dict
from services.filter_bam_multifile import create_read_dict
from services.filter_bam_multifile import filter_reads
from services.filter_bam_multifile import create_dict_of_transcripts_and_reads
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

    def generate_reads_and_locations(self, dict_of_transcripts_and_reads: dict):
        reads_and_locations = {}
        for key, value in self.matching_cases_dict.items():
            if not key[0] in dict_of_transcripts_and_reads:
                continue
            for read in dict_of_transcripts_and_reads[key[0]]:
                if read not in reads_and_locations:
                    reads_and_locations[read] = []
                reads_and_locations[read].append(value)
        return reads_and_locations

    def execute(self):
        output_manager.output_line({
            "line": "PROCESSING BAM-FILE",
            "is_title": True
        })
        output_manager.output_line({
            "line": f"Input BAM-file: {self.bam_path}",
            "is_info": True
        })

        # output_filename_dict = create_output_filename_dict(
        #     self.bam_path,
        #     self.transcript_set,
        #     TEMPORARY_DIR
        # )
        # read_dict = create_read_dict(
        #     output_filename_dict,
        #     self.tsv_path
        # )
        # filter_reads(
        #     self.bam_path,
        #     output_filename_dict,
        #     read_dict
        # )

        # self.iterate_extracted_files()

        output_manager.output_line({
            "line": "Extracting reads from tsv-file",
            "is_info": True
        })

        dict_of_transcripts_and_reads = create_dict_of_transcripts_and_reads(
            self.transcript_set, self.tsv_path)

        reads_and_locations = self.generate_reads_and_locations(
            dict_of_transcripts_and_reads)

        output_manager.output_line({
            "line": "NUMBER OF MATCHING CASES:" + str(len(self.matching_cases_dict)),
            "is_info": True
        })

        output_manager.output_line({
            "line": "NUMBER OF READS: " + str(len(reads_and_locations)),
            "is_info": True
        })
        alignment_parser.execute(self.bam_path, reads_and_locations)
        output_manager.output_line({
            "line": "Insertions and deletions found at given locations",
            "is_info": True
        })
        for line in alignment_parser.case_count:
            output_manager.output_line({
                "line": line,
                "is_info": True
            })

        self.remove_temporary_path()

    def create_temporary_path(self):
        if not os.path.exists(TEMPORARY_DIR):
            os.mkdir(TEMPORARY_DIR)

    def remove_temporary_path(self):
        if os.path.exists(TEMPORARY_DIR):
            for file in os.listdir(TEMPORARY_DIR):
                os.remove(os.path.join(TEMPORARY_DIR, file))
            os.rmdir(TEMPORARY_DIR)

    def iterate_extracted_files(self):
        files = {
            "found": 0,
            "not_found": 0
        }
        for key, value in self.matching_cases_dict.items():
            filename = Path(self.bam_path).stem + "." + key[0] + ".bam"

            if filename in os.listdir(TEMPORARY_DIR):
                files["found"] += 1
            else:
                files["not_found"] += 1
            if filename in os.listdir(TEMPORARY_DIR):
                alignment_parser.execute(filename, location=value)

        print(files)
        output_manager.output_line({
            "line": "Insertions and deletions found at given locations",
            "is_info": True
        })

        for key, value in alignment_parser.case_count.items():
            output_manager.output_line({
                "line": f"{key}: {value}",
                "is_info": True
            })
