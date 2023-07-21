from pyfaidx import Fasta
from services.output_manager import default_output_manager as output_manager
from config import DEFAULT_WINDOW_SIZE


class FastaExtractor:

    def __init__(self, fasta_config: dict):
        self.fasta = None
        self.fasta_path = fasta_config.get('fasta_path', '')
        self.offset = fasta_config.get('offset', -1)
        self.matching_cases_dict = fasta_config.get(
            'matching_cases_dict', {})
        self.window_size = int(fasta_config.get(
            'window_size', DEFAULT_WINDOW_SIZE))
        self.closest_canonicals_overview = {"left": {}, "right": {}}
        self.index_correction = -1

    def initialize_fasta(self):
        if self.fasta_path:
            self.fasta = Fasta(self.fasta_path)

    def extract_characters_at_given_coordinates(self, coordinates: tuple):
        if self.fasta:
            chromosome, start, end = coordinates
            return self.fasta[chromosome][start + self.index_correction:end + self.index_correction]
        return -1

    def output_section_header(self):
        output_manager.output_line({
            "line": "CLOSEST CANONICALS",
            "is_title": True
        })
        output_manager.output_line({
            "line": "Fethching reference FASTA-file.",
            "is_info": True
        })

    def check_errors(self) -> bool:
        errors = False
        if not self.fasta:
            output_manager.output_line({
                "line": "FASTA-file not found. Please check path and try again. " +
                "Moving to next section.",
                "is_error": True
            })
            errors = True

        if self.offset == -1:
            output_manager.output_line({
                "line": "No offset value given. Nothing to do here.",
                "is_error": True
            })
            errors = True
        if not self.matching_cases_dict:
            output_manager.output_line({
                "line": "No matching cases found. Nothing to do here.",
                "is_error": True
            })
        return errors

    def find_closest_canonicals(self, nucleotides: str, dict_key: tuple, canonicals: list):
        nucleotides_middle = int(len(nucleotides) / 2)
        closest_canonicals = {}
        aligned_splice_site_nucleotides = nucleotides[nucleotides_middle:nucleotides_middle + 2]
        for i in range(1, nucleotides_middle):
            left_position = nucleotides[nucleotides_middle -
                                        i:nucleotides_middle - i + 2]
            right_position = nucleotides[nucleotides_middle +
                                         i:nucleotides_middle + i + 2]
            if left_position in canonicals and 'left' not in closest_canonicals:
                closest_canonicals['left'] = (
                    left_position, aligned_splice_site_nucleotides, i)
            if right_position in canonicals and 'right' not in closest_canonicals:
                closest_canonicals['right'] = (
                    right_position, aligned_splice_site_nucleotides, i)
        for item in ['left', 'right']:
            if item not in closest_canonicals:
                closest_canonicals[item] = (
                    aligned_splice_site_nucleotides, aligned_splice_site_nucleotides, 0)
        self.matching_cases_dict[dict_key]["closest_canonical"] = closest_canonicals

        left_result = self.matching_cases_dict[dict_key]["closest_canonical"]["left"]
        right_result = self.matching_cases_dict[dict_key]["closest_canonical"]["right"]
        if left_result not in self.closest_canonicals_overview["left"]:
            self.closest_canonicals_overview["left"][left_result] = 0
        self.closest_canonicals_overview["left"][left_result] += 1
        if right_result not in self.closest_canonicals_overview["right"]:
            self.closest_canonicals_overview["right"][right_result] = 0
        self.closest_canonicals_overview["right"][right_result] += 1

    def iterate_matching_cases(self):
        for key, value in self.matching_cases_dict.items():
            if value["location_type"] == "start":
                splice_cite_location = value["location"] - 2
            else:
                splice_cite_location = value["location"] + 1
            coordinates = (value['seq_id'],
                           splice_cite_location - self.window_size,
                           splice_cite_location + self.window_size)
            nucleotides = self.extract_characters_at_given_coordinates(
                coordinates)
            if value["location_type"] == "start":
                canonicals = ["AG", "AC"]
            else:
                canonicals = ["GT", "GC", "AT"]
            self.find_closest_canonicals(str(nucleotides), key, canonicals)

    def execute_fasta_extraction(self):
        self.output_section_header()
        self.initialize_fasta()
        if self.check_errors():
            return

        output_manager.output_line({
            "line": "Searching closest possible canonical " +
            f"sites for provided offset range: {self.offset}",
            "is_info": True
        })

        self.iterate_matching_cases()

        output_manager.output_line({
            "line": "Finished.",
            "is_info": True
        })
