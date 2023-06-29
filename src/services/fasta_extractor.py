import json
from pyfaidx import Fasta
from services.output_manager import default_output_manager as output_manager
from config import FASTA_OVERVIEW_FILE, DEFAULT_WINDOW_SIZE


class FastaExtractor:

    def __init__(self, fasta_config: dict):
        self.fasta = None
        self.fasta_path = fasta_config.get('fasta_path', '')
        self.offset = fasta_config.get('offset', -1)
        self.gffcompare_gtf = fasta_config.get('gffcompare_gtf', '')
        self.reference_gtf = fasta_config.get('reference_gtf', '')
        self.class_codes = fasta_config.get('class_codes', '')
        self.matching_cases_dict = fasta_config.get(
            'matching_cases_dict', {})
        self.window_size = int(fasta_config.get(
            'window_size', DEFAULT_WINDOW_SIZE))
        self.closest_canonicals_overview = {"left": {}, "right": {}}

    def initialize_fasta(self):
        if self.fasta_path:
            self.fasta = Fasta(self.fasta_path)

    def extract_characters_at_given_coordinates(self, coordinates: tuple):
        if self.fasta:
            chromosome, start, end = coordinates
            return self.fasta[chromosome][start:end]
        return -1

    def output_section_header(self):
        output_manager.output_line({
            "line": "FASTA EXTRACTION",
            "is_title": True
        })
        output_manager.output_line({
            "line": "Fethching reference fasta file...",
            "is_info": True
        })

    def check_errors(self) -> bool:
        errors = False
        if not self.fasta:
            output_manager.output_line({
                "line": "fasta-file not found. Please check path and try again. Moving to next section.",
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

    def write_results_to_file(self):
        json_overview = {
            "strand": {
                "+": {
                    "start": {},
                    "end": {}
                },
                "-": {
                    "start": {},
                    "end": {}
                }
            },

        }
        for key, value in self.matching_cases_dict.items():
            strand, location_type = value['strand'], value['location_type']
            splice_cite_sequence = str(value['splice_cite_sequence'])
            if splice_cite_sequence not in json_overview['strand'][strand][location_type]:
                json_overview['strand'][strand][location_type][splice_cite_sequence] = 0
            json_overview['strand'][strand][location_type][splice_cite_sequence] += 1

        for strand in json_overview['strand']:
            for position in json_overview['strand'][strand]:
                json_overview['strand'][strand][position] = dict(
                    sorted(
                        json_overview['strand'][strand][position].items(),
                        key=lambda item: item[1],
                        reverse=True
                    )
                )

        with open(FASTA_OVERVIEW_FILE, "w", encoding="utf-8") as file:
            file.write("# Overview\n")
            file.write("## Offset characters\n")
            file.write(
                "This section contains the characters in the reference ")
            file.write(
                "FASTA-file at the offset specified by the user.  \n\n")
            file.write("**Arguments provided by the user:**\n")
            file.write("```\n")
            file.write("gffcompare GTF-file:\n" +
                       self.gffcompare_gtf + "\n\n")
            file.write("Reference GTF-file:\n" +
                       self.reference_gtf + "\n\n")
            file.write("Reference FASTA-file:\n" +
                       self.fasta_path + "\n\n")
            file.write("Specified offset: " + str(self.offset) + "\n")
            file.write("Class codes: " + str(self.class_codes) + "\n")
            file.write("```\n")
            file.write("**Results in JSON-format:**  \n")
            file.write("```json\n")
            file.write(json.dumps(json_overview, indent=4))
            file.write("\n```\n")
            file.write("\n\n## Results in table-format \n")
            file.write("This section contains the results in table-format.  \n")
            current_chromosome = "chr1"
            file.write(f'\n### {current_chromosome}\n')
            file.write(
                "| transcript | read_id | strand | exon | nucleotides |\n")
            file.write("| --- | --- | --- | --- | --- |\n")
            for key, value in self.matching_cases_dict.items():
                chromosome = value['transcript_id'].split('.')[1]
                if current_chromosome != chromosome:
                    file.write(f'\n### {chromosome}\n')
                    file.write(
                        "| transcript | strand | exon | nucleotides |\n")
                    file.write("| --- | --- | --- | --- | --- |\n")
                    current_chromosome = chromosome
                file.write("| " + str(value['transcript_id']) + " | " +
                           " | " + str(value['strand']) + " | " + str(value['exon_number']) +
                           " | " + str(value['splice_cite_sequence']) + " |\n")
            file.close()
            output_manager.output_line({
                "line": "Results written to file: " + FASTA_OVERVIEW_FILE,
                "is_info": True
            })

    def extract_splice_cite_sequence(self):
        for key, value in self.matching_cases_dict.items():
            chromosome = value['transcript_id'].split('.')[1]
            location = value['location']
            if value['location_type'] == 'end':
                coordinates = (chromosome, location - 1, location + 1)
            else:
                coordinates = (chromosome, location, location + 2)
            chars = self.extract_characters_at_given_coordinates(
                coordinates)
            value['splice_cite_sequence'] = str(chars)

        self.write_results_to_file()

    def extract_window_nucleotides(self):

        for key, value in self.matching_cases_dict.items():
            chromosome = value['transcript_id'].split('.')[1]
            location = value['location']
            if value['location_type'] == 'end':
                coordinates = (chromosome, location -
                               self.window_size - 1, location + 1)
            else:
                coordinates = (chromosome, location,
                               location + self.window_size)
            nucleotides = self.extract_characters_at_given_coordinates(
                coordinates)
            self.matching_cases_dict[key]['position_sequence'] = str(
                nucleotides)

    def find_closest_canonicals(self, nucleotides: str, dict_key: str, canonicals: list):
        nucleotides_middle = int(len(nucleotides) / 2)
        closest_canonicals = {}
        aligned_splice_site_nucleotides = nucleotides[nucleotides_middle:nucleotides_middle + 2]
        for i in range(1, nucleotides_middle):
            if nucleotides[nucleotides_middle - i:nucleotides_middle - i + 2] in canonicals and 'left' not in closest_canonicals:
                closest_canonicals['left'] = (
                    nucleotides[nucleotides_middle - i:nucleotides_middle - i + 2], aligned_splice_site_nucleotides)
            if nucleotides[nucleotides_middle + i:nucleotides_middle + i + 2] in canonicals and 'right' not in closest_canonicals:
                closest_canonicals['right'] = (
                    nucleotides[nucleotides_middle + i:nucleotides_middle + i + 2], aligned_splice_site_nucleotides)
        if 'left' not in closest_canonicals:
            closest_canonicals['left'] = (
                aligned_splice_site_nucleotides, aligned_splice_site_nucleotides)
        if 'right' not in closest_canonicals:
            closest_canonicals['right'] = (
                aligned_splice_site_nucleotides, aligned_splice_site_nucleotides)
        self.matching_cases_dict[dict_key]["closest_canonicals"] = closest_canonicals

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
            coordinates = (key.split('.')[
                           1], value["location"] - self.window_size, value["location"] + self.window_size)
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
            "line": f"Closed offset range: {self.offset}",
            "is_info": True
        })

        self.extract_splice_cite_sequence()
        self.extract_window_nucleotides()
