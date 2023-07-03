import os

from services.read_management import create_dict_of_transcripts_and_reads
from services.output_manager import default_output_manager as output_manager
from services.alignment_parser import default_alignment_parser as alignment_parser
from services.log_manager import default_log_manager as log_manager

from config import LOG_FILE_DIR


class BamManager:

    def __init__(self,
                 bam_path: str,
                 tsv_path: str,
                 matching_cases_dict: dict,
                 extended_debugging: bool = False):
        self.bam_path = bam_path
        self.tsv_path = tsv_path
        self.matching_cases_dict = matching_cases_dict
        self.transcript_set = set()
        for value in self.matching_cases_dict.values():
            self.transcript_set.add(value['transcript_id'])
        self.extended_debugging = extended_debugging
        self.debug_log_path_transcripts_and_reads = os.path.join(
            LOG_FILE_DIR, "dict_of_transcripts_and_reads.log")
        self.debug_log_path_reads_and_locations = os.path.join(
            LOG_FILE_DIR, "reads_and_locations.log")

    def write_debug_logs(self, transcripts_and_reads: dict, reads_and_locations: dict):

        with open(self.debug_log_path_transcripts_and_reads, "w") as f:
            for key, value in transcripts_and_reads.items():
                f.write(f"{key}\t{value}\n")

        with open(self.debug_log_path_reads_and_locations, "w") as f:
            for key, value in reads_and_locations.items():
                f.write(f"{key}\t{value}\n")

    def generate_reads_and_references(self, dict_of_transcripts_and_reads: dict):
        reads_and_references = {}
        for matching_case_key, matching_case_values in self.matching_cases_dict.items():
            if not matching_case_values['transcript_id'] in dict_of_transcripts_and_reads:
                continue
            for read in dict_of_transcripts_and_reads[matching_case_values['transcript_id']]:
                if read not in reads_and_references:
                    reads_and_references[read] = set()
                reads_and_references[read].add(matching_case_key)
        return reads_and_references

    def output_heading_information(self):
        output_manager.output_line({
            "line": "PROCESSING BAM-FILE",
            "is_title": True
        })
        output_manager.output_line({
            "line": f"Input BAM-file: {self.bam_path}",
            "is_info": True
        })

        output_manager.output_line({
            "line": "Extracting reads from tsv-file",
            "is_info": True
        })

    def generate_dictionaries(self):
        dict_of_transcripts_and_reads = create_dict_of_transcripts_and_reads(
            self.transcript_set, self.tsv_path)

        reads_and_references = self.generate_reads_and_references(
            dict_of_transcripts_and_reads)

        log_manager.debug_logs['transcripts_and_reads'] = dict_of_transcripts_and_reads
        log_manager.debug_logs['reads_and_references'] = reads_and_references

        output_manager.output_line({
            "line": "NUMBER OF MATCHING CASES:" + str(len(self.matching_cases_dict)),
            "is_info": True
        })

        output_manager.output_line({
            "line": "NUMBER OF READS: " + str(len(reads_and_references)),
            "is_info": True
        })

        output_manager.output_line({
            "line": "Analyzing offset of reads. This may take a while.",
            "is_info": True
        })

        return reads_and_references, dict_of_transcripts_and_reads

    def execute(self, window_size: int):

        self.output_heading_information()

        reads_and_references, dict_of_transcripts_and_reads = self.generate_dictionaries()

        alignment_parser.execute(
            self.bam_path,
            window_size,
            reads_and_references,
            self.matching_cases_dict,
            dict_of_transcripts_and_reads
        )
