import os
import pysam
from services.output_manager import default_output_manager as output_manager
from config import DEFAULT_WINDOW_SIZE, LOG_FILE_DIR


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
            Follows the singleton pattern. Total case count for all processed BAM-files
            is stored in self.case_count.
        """

        self.reads_and_transcripts = {}
        self.case_count = {
            "insertions": {'+': {}, '-': {}},
            "deletions": {'+': {}, '-': {}},
        }
        self.updated_case_count = {}
        # TODO: add to configuration or to user arguments (window size)
        self.window_size = int(DEFAULT_WINDOW_SIZE)
        self.error_file_output_dir = os.path.join(
            LOG_FILE_DIR, "alignment_errors.log")

    def initialize_file(self, filename: str):
        # pylint: disable=no-member
        self.samfile = pysam.AlignmentFile(filename, "rb")

    def extract_location_from_cigar_string(self, cigar_tuples: list, reference_start: int, reference_end: int, location: int):
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
        deletions = 0
        insertions = 0
        debug_list = []
        location = 0

        cigar_code_list = []

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

        del_key = ("deletions", strand, loc_type, offset)
        if del_key not in self.updated_case_count:
            self.updated_case_count[del_key] = {}
        if deletions not in self.updated_case_count[del_key]:
            self.updated_case_count[del_key][deletions] = 0
        self.updated_case_count[del_key][deletions] += 1
        if deletions not in self.case_count["deletions"][strand]:
            self.case_count["deletions"][strand][deletions] = 0
        self.case_count["deletions"][strand][deletions] += 1

        ins_key = ("insertions", strand, loc_type, offset)
        if ins_key not in self.updated_case_count:
            self.updated_case_count[ins_key] = {}
        if insertions not in self.updated_case_count[ins_key]:
            self.updated_case_count[ins_key][insertions] = 0
        self.updated_case_count[ins_key][insertions] += 1
        if insertions not in self.case_count["insertions"][strand]:
            self.case_count["insertions"][strand][insertions] = 0
        self.case_count["insertions"][strand][insertions] += 1

        if deletions >= self.window_size or insertions >= self.window_size:
            debug_list.append(
                f"deletions: {deletions}, insertions: {insertions}")
            return True, debug_list
        return False, debug_list

    def write_alignment_errors_to_file(self, errors: list):
        with open(self.error_file_output_dir, "w") as file:
            file.write(
                "qname\ttranscripts\tlocation\talign_location\ttype\tread.reference_start\tread.reference_end\tlist of alignments\n")
            file.writelines(errors)

    def process_bam_file(self, reads_and_locations: dict):
        count = 0
        errors = []
        set_of_processed_reads = set()
        for read in self.samfile.fetch():
            if read.query_name not in set_of_processed_reads:
                set_of_processed_reads.add(read.query_name)
            else:
                output_manager.output_line({
                    "line": "Error: read " + str(read.query_name) + "already processed",
                    "is_error": True
                })
                continue

            if read.is_supplementary:
                continue
            if read.query_name in reads_and_locations:
                count += 1
                if count % 1000 == 0:
                    output_manager.output_line({
                        "line": "Processed " + str(count) + " reads",
                        "end_line": "\r",
                        "is_info": True
                    })
                for dict_item in reads_and_locations[read.query_name]:
                    location, loc_type = dict_item["location"], dict_item["location_type"]
                    strand = dict_item["strand"]
                    offset = dict_item["offset"]
                    idx_corrected_location = location - 1

                    if read.reference_start > idx_corrected_location or read.reference_end < idx_corrected_location:
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

                    response, debug_list = self.count_indels_from_cigar_codes_in_given_window(
                        read.cigartuples,
                        aligned_location,
                        loc_type,
                        strand,
                        offset)

                    if response:
                        errors.append(
                            f"{read.query_name}\t{self.reads_and_transcripts[read.query_name]}\t{idx_corrected_location}\t{aligned_location}\t{loc_type}\t{read.reference_start}\t{read.reference_end}\t{debug_list}\n")
        if errors:
            self.write_alignment_errors_to_file(errors)

        output_manager.output_line({
            "line": "\nFinished",
            "is_info": True
        })

    # TODO: remove transcripts and reads. It is here for debugging
    def execute(self, filename: str, window_size: int,
                reads_and_locations: dict, transcripts_and_reads: dict):
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
        self.process_bam_file(reads_and_locations)


default_alignment_parser = AlignmentParser()
