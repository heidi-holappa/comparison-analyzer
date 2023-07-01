import pysam

from services.output_manager import default_output_manager as output_manager
from services.log_manager import default_log_manager as log_manager
from config import DEFAULT_WINDOW_SIZE


class AlignmentParser:

    def __new__(cls):
        """
            Singleton pattern.
        """
        if not hasattr(cls, 'instance'):
            cls.instance = super(AlignmentParser, cls).__new__(cls)
        return cls.instance

    # TODO: remove reads and transcripts. It is for debugging purposes.
    def __init__(self):
        """
            Parse a BAM-file and count the number of insertions and deletions at a given location.
            Follows the singleton pattern. 
        """

        self.reads_and_transcripts = {}
        # TODO: add to configuration or to user arguments (window size)
        self.window_size = int(DEFAULT_WINDOW_SIZE)

    def initialize_file(self, filename: str):
        # pylint: disable=no-member
        self.samfile = pysam.AlignmentFile(filename, "rb")

    def extract_location_from_cigar_string(self,
                                           cigar_tuples: list,
                                           reference_start: int,
                                           reference_end: int,
                                           location: int):
        relative_position = location - reference_start
        alignment_position = 0
        ref_position = 0

        for idx, cigar_code in enumerate(cigar_tuples):

            if cigar_code[0] in [0, 2, 3, 7, 8]:
                ref_position += cigar_code[1]
            if ref_position <= relative_position and not reference_start + ref_position == reference_end:
                alignment_position += cigar_code[1]
            else:
                return alignment_position + (cigar_code[1] - (ref_position - relative_position))

        return -1

    def count_indels_from_cigar_codes_in_given_window(self,
                                                      cigar_tuples: list,
                                                      aligned_location: int,
                                                      loc_type: str,
                                                      strand: str,
                                                      offset: int = 0):
        """
        Get cigar codes in a given window.

        Args:
            cigar_tuples (list): list of cigar tuples (cigar code, aligned position)
            aligned_location (int): aligned location
            loc_type (str): type of location (start or end)
        """
        result = {
            'deletions': 0,
            'insertions': 0
        }

        deletions = 0
        insertions = 0

        cigar_code_list = []
        location = 0

        debug_list = []

        if loc_type == "end":
            aligned_location = aligned_location - self.window_size + 1

        for cigar_code in cigar_tuples:
            if self.window_size == len(cigar_code_list):
                break
            if location + cigar_code[1] > aligned_location:
                overlap = location + \
                    cigar_code[1] - (aligned_location + len(cigar_code_list))
                cigar_code_list.extend(
                    [cigar_code[0] for _ in range(min(self.window_size - len(cigar_code_list), overlap))])
            location += cigar_code[1]

        debug_list.append(cigar_code_list)
        for cigar_code in cigar_code_list:
            if cigar_code == 2:
                deletions += 1
            if cigar_code == 1:
                insertions += 1

        result['deletions'] = deletions
        result['insertions'] = insertions

        if deletions >= self.window_size or insertions >= self.window_size:
            debug_list.append(
                f"deletions: {deletions}, insertions: {insertions}")
            return True, debug_list, result
        return False, debug_list, result

    def process_bam_file(self,
                         reads_and_references: dict,
                         matching_cases_dict: dict):
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
                    location = matching_cases_dict[matching_case_key]["location"]
                    loc_type = matching_cases_dict[matching_case_key]["location_type"]
                    strand = matching_cases_dict[matching_case_key]["strand"]
                    offset = matching_cases_dict[matching_case_key]["offset"]

                    idx_corrected_location = location - 1

                    if read.reference_start > idx_corrected_location or read.reference_end < idx_corrected_location:
                        errors.append(
                            f"Non-matching location: {read.query_name}, {matching_case_key}\t")
                        continue

                    if not read.cigartuples:
                        continue

                    if not read.reference_end:
                        errors.append(f"{read.query_name}\tno reference end\n")
                        continue

                    aligned_location = self.extract_location_from_cigar_string(
                        read.cigartuples,
                        read.reference_start,
                        read.reference_end,
                        idx_corrected_location
                    )

                    error, debug_list, result = self.count_indels_from_cigar_codes_in_given_window(
                        read.cigartuples,
                        aligned_location,
                        loc_type,
                        strand,
                        offset)
                    matching_cases_dict[matching_case_key]['indel_errors'] = result

                    if error:
                        log_manager.alignment_erros.append(
                            f"{read.query_name}\t{self.reads_and_transcripts[read.query_name]}\t" +
                            f"{idx_corrected_location}\t{aligned_location}\t{loc_type}\t" +
                            f"{read.reference_start}\t{read.reference_end}\t{debug_list}\n")

        if prev_processed_reads_counter:
            output_manager.output_line({
                "line": str(prev_processed_reads_counter) +
                " iterations extracted an already processed read from the BAM-file",
                "is_error": True,
            })

        output_manager.output_line({
            "line": "\nFinished",
            "is_info": True
        })

    # TODO: remove transcripts and reads. It is here for debugging
    def execute(self,
                filename: str,
                window_size: int,
                reads_and_references: dict,
                matching_cases_dict: dict,
                transcripts_and_reads: dict):
        """
        Initialize AlignmentFile and execute the alignment parser.

        Args:
            filename (str): path of BAM-file to iterate
            location (int): coordinate of the interesting event
        """

        self.window_size = int(window_size)

        for key, value in transcripts_and_reads.items():
            for read in value:
                if read not in self.reads_and_transcripts:
                    self.reads_and_transcripts[read] = set()
                self.reads_and_transcripts[read].add(key)

        self.initialize_file(filename)
        self.process_bam_file(reads_and_references, matching_cases_dict)


default_alignment_parser = AlignmentParser()
