import json
from pyfaidx import Fasta
from services.output_manager import default_output_manager as output_manager


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

    def initialize_fasta(self):
        if self.fasta_path:
            self.fasta = Fasta(self.fasta_path)

    def extract_characters_at_given_coordinates(self, coordinates: tuple):
        if self.fasta:
            chromosome, start, end = coordinates
            return self.fasta[chromosome][start:end]

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

    def write_results_to_file(self, results: dict):
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
        for key, value in results.items():
            if str(value) not in json_overview['strand'][key[2]][key[4]]:
                json_overview['strand'][key[2]][key[4]][str(value)] = 0
            json_overview['strand'][key[2]][key[4]][str(value)] += 1

        for strand in json_overview['strand']:
            for position in json_overview['strand'][strand]:
                json_overview['strand'][strand][position] = dict(
                    sorted(
                        json_overview['strand'][strand][position].items(),
                        key=lambda item: item[1],
                        reverse=True
                    )
                )

        with open("overview.md", "w", encoding="utf-8") as file:
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

    def execute_fasta_extraction(self):
        self.output_section_header()
        self.initialize_fasta()
        if self.check_errors():
            return

        output_manager.output_line({
            "line": f"Offset value: {self.offset}",
            "is_info": True
        })

        results = {}
        for key, value in self.matching_cases_dict.items():
            chromosome = key[0].split('.')[1]
            if key[4] == 'end':
                coordinates = (chromosome, value, value + 2)
            else:
                coordinates = (chromosome, value - 3, value - 1)
            chars = self.extract_characters_at_given_coordinates(
                coordinates)
            results[key] = chars

        self.write_results_to_file(results)
        return results
