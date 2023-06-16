import os

from services.read_management import create_dict_of_transcripts_and_reads
from services.output_manager import default_output_manager as output_manager
from services.alignment_parser import default_alignment_parser as alignment_parser
from services.graph_manager import default_graph_manager as graph_manager

from config import TEMPORARY_DIR, LOG_FILE_DIR


class BamManager:

    def __init__(self, bam_path: str, tsv_path: str, matching_cases_dict: dict):
        self.bam_path = bam_path
        self.tsv_path = tsv_path
        self.matching_cases_dict = matching_cases_dict
        self.transcript_set = set()
        for row in matching_cases_dict:
            self.transcript_set.add(row[0])

    def write_debug_logs(self, transcripts_and_reads: dict, reads_and_locations: dict):
        dict_of_transcripts_and_reads_log = os.path.join(
            LOG_FILE_DIR, "dict_of_transcripts_and_reads.log")

        with open(dict_of_transcripts_and_reads_log, "w") as f:
            for key, value in transcripts_and_reads.items():
                f.write(f"{key}\t{value}\n")

        reads_and_locations_log = os.path.join(
            LOG_FILE_DIR, "reads_and_locations.log")
        with open(reads_and_locations_log, "w") as f:
            for key, value in reads_and_locations.items():
                f.write(f"{key}\t{value}\n")

    def generate_reads_and_locations(self, dict_of_transcripts_and_reads: dict):
        reads_and_locations = {}
        for key, value in self.matching_cases_dict.items():
            if not key[0] in dict_of_transcripts_and_reads:
                continue
            for read in dict_of_transcripts_and_reads[key[0]]:
                if read not in reads_and_locations:
                    reads_and_locations[read] = []
                location_and_type = (value, key[4])
                reads_and_locations[read].append(location_and_type)
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

        output_manager.output_line({
            "line": "Extracting reads from tsv-file",
            "is_info": True
        })

        dict_of_transcripts_and_reads = create_dict_of_transcripts_and_reads(
            self.transcript_set, self.tsv_path)

        reads_and_locations = self.generate_reads_and_locations(
            dict_of_transcripts_and_reads)

        self.write_debug_logs(
            dict_of_transcripts_and_reads, reads_and_locations)

        output_manager.output_line({
            "line": "NUMBER OF MATCHING CASES:" + str(len(self.matching_cases_dict)),
            "is_info": True
        })

        output_manager.output_line({
            "line": "NUMBER OF READS: " + str(len(reads_and_locations)),
            "is_info": True
        })

        output_manager.output_line({
            "line": "Analyzing offset of reads. This may take a while.",
            "is_info": True
        })

        alignment_parser.execute(
            self.bam_path, reads_and_locations, dict_of_transcripts_and_reads)

        output_manager.output_line({
            "line": "Insertions and deletions found at given locations",
            "is_info": True
        })
        print(alignment_parser.case_count)
        for key, value in alignment_parser.case_count.items():
            output_manager.output_line({
                "line": f"{key}: {sorted(value.items())}",
                "is_info": True
            })
            if value:
                graph_manager.construct_bar_chart_from_dict(
                    graph_values=value,
                    title=key,
                    x_label="Number of cases",
                    y_label="Number of reads",
                )

        self.remove_temporary_path()

    def create_temporary_path(self):
        if not os.path.exists(TEMPORARY_DIR):
            os.mkdir(TEMPORARY_DIR)

    def remove_temporary_path(self):
        if os.path.exists(TEMPORARY_DIR):
            for file in os.listdir(TEMPORARY_DIR):
                os.remove(os.path.join(TEMPORARY_DIR, file))
            os.rmdir(TEMPORARY_DIR)
