import pysam

from services.read_management import create_dict_of_transcripts_and_reads
from services.output_manager import default_output_manager as output_manager
from services.alignment_parser import default_alignment_parser as alignment_parser
from services.log_manager import default_log_manager as log_manager


class BamManager:

    def __init__(self,
                 bam_path: str,
                 tsv_path: str,
                 matching_cases_dict: dict):
        # pylint: disable=no-member

        self.bam_path = bam_path
        self.tsv_path = tsv_path
        self.samfile = pysam.AlignmentFile(self.bam_path, "rb")

        self.matching_cases_dict = matching_cases_dict
        self.reads_and_transcripts = {}

    def process_bam_file(self, reads_and_references: dict, window_size: int):
        output_manager.output_line({
            "line": "Iterating reads and counting indels. This may take a while.",
            "is_info": True
        })
        count = 0
        errors = []
        set_of_processed_reads = set()
        prev_processed_reads_counter = 0
        for read in self.samfile.fetch():

            if read.is_supplementary:
                continue
            if read.query_name in reads_and_references:
                if read.query_name not in set_of_processed_reads:
                    set_of_processed_reads.add(read.query_name)
                else:
                    prev_processed_reads_counter += 1
                    continue
                count += 1
                if count % 1000 == 0:
                    output_manager.output_line({
                        "line": "Processed " + str(count) + " reads",
                        "end_line": "\r",
                        "is_info": True,
                        "save_to_log": False
                    })
                for matching_case_key in reads_and_references[read.query_name]:
                    location = self.matching_cases_dict[matching_case_key]["location"]
                    loc_type = self.matching_cases_dict[matching_case_key]["location_type"]

                    if 'indel_errors' not in self.matching_cases_dict[matching_case_key]:
                        self.matching_cases_dict[matching_case_key]['indel_errors'] = {
                            'insertions': {},
                            'deletions': {}
                        }

                    idx_corrected_location = location - 1

                    if read.reference_start > idx_corrected_location \
                            or read.reference_end < idx_corrected_location:
                        errors.append(
                            f"Non-matching location: {read.query_name}, {matching_case_key}\t")
                        continue

                    if not read.cigartuples:
                        continue

                    if not read.reference_end:
                        errors.append(f"{read.query_name}\tno reference end\n")
                        continue

                    aligned_location = alignment_parser.extract_location_from_cigar_string(
                        read.cigartuples,
                        read.reference_start,
                        read.reference_end,
                        idx_corrected_location
                    )

                    error, debug_list, result = alignment_parser.count_indels_from_cigar_codes_in_given_window(
                        read.cigartuples,
                        aligned_location,
                        loc_type,
                        window_size)
                    for key, value in result.items():
                        if value not in self.matching_cases_dict[matching_case_key]['indel_errors'][key]:
                            self.matching_cases_dict[matching_case_key]['indel_errors'][key][value] = 0
                        self.matching_cases_dict[matching_case_key]['indel_errors'][key][value] += 1

                    if error:
                        log_manager.alignment_erros.append(
                            f"{read.query_name}\t{self.reads_and_transcripts[read.query_name]}\t" +
                            f"{idx_corrected_location}\t{aligned_location}\t{loc_type}\t" +
                            f"{read.reference_start}\t{read.reference_end}\t{debug_list}\n")

        if prev_processed_reads_counter:
            output_manager.output_line({
                "line": str(prev_processed_reads_counter) +
                " iterations extracted an already processed read from the BAM-file",
                "is_warning": True,
            })

        output_manager.output_line({
            "line": "Processing BAM-file finished.",
            "is_info": True
        })

    def output_heading_information(self):
        output_manager.output_line({
            "line": "ITERATE READS",
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

    def generate_reads_and_transcripts(self, transcripts_and_reads: dict):
        for key, value in transcripts_and_reads.items():
            for read in value:
                if read not in self.reads_and_transcripts:
                    self.reads_and_transcripts[read] = set()
                self.reads_and_transcripts[read].add(key)

    def generate_dictionaries(self):
        transcript_set = set()
        for value in self.matching_cases_dict.values():
            transcript_set.add(value['transcript_id'])

        transcripts_and_reads = create_dict_of_transcripts_and_reads(
            transcript_set, self.tsv_path)

        reads_and_references = self.generate_reads_and_references(
            transcripts_and_reads)

        self.generate_reads_and_transcripts(transcripts_and_reads)

        log_manager.debug_logs['transcripts_and_reads'] = transcripts_and_reads
        log_manager.debug_logs['reads_and_references'] = reads_and_references

        cases_line = f"Number of matching cases: {len(self.matching_cases_dict)}, " + \
            f"number of reads: {len(reads_and_references)}\n"

        output_manager.output_line({
            "line": cases_line,
            "is_info": True
        })

        output_manager.output_line({
            "line": "Note: one read may be related to multiple matching cases.",
            "is_info": True
        })

        return reads_and_references

    def execute(self, window_size: int):

        self.output_heading_information()

        reads_and_references = self.generate_dictionaries()

        self.process_bam_file(reads_and_references, window_size)
