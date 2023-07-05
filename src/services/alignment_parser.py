import pysam

from services.output_manager import default_output_manager as output_manager
from services.log_manager import default_log_manager as log_manager
from config import DEFAULT_WINDOW_SIZE


class AlignmentParser:

    def __init__(self):
        """
            Parse a BAM-file and count the number of insertions and deletions at a given location.
            Follows the singleton pattern. 
        """

        # TODO: add to configuration or to user arguments (window size)
        self.window_size = int(DEFAULT_WINDOW_SIZE)

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
                                                      loc_type: str):
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


default_alignment_parser = AlignmentParser()
