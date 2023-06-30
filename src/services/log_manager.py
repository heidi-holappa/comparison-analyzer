import json
from services.output_manager import default_output_manager as output_manager
from services.graph_manager import default_graph_manager as graph_manager
from config import LOG_FILE_DIR, FASTA_OVERVIEW_FILE


class LogManager:

    def __init__(self, matching_cases_dict: dict, parser_args):
        self.matching_cases_dict = matching_cases_dict
        self.parser_args = parser_args

        pass

    def compute_json_overview_dict_for_closest_canonicals(self):
        json_overview = {}
        for value in self.matching_cases_dict.values():
            strand, location_type = value['strand'], value['location_type']
            offset = value['offset']
            for indel_errors_key in value['indel_errors']:
                json_key = (indel_errors_key, strand, location_type, offset)
                for canonicals_key, canonicals_value in value['closest_canonical'].items():
                    if json_key not in json_overview:
                        json_overview[json_key] = {}
                    if canonicals_key not in json_overview[json_key]:
                        json_overview[json_key][canonicals_key] = {}
                    for canonical_pair in canonicals_key:
                        if canonical_pair not in json_overview[json_key][canonicals_key]:
                            json_overview[json_key][canonicals_key][canonical_pair] = 0
                        json_overview[json_key][canonicals_key][canonical_pair] += 1

        for site_locations in json_overview.values():
            for results in site_locations.values():
                results = dict(
                    sorted(
                        results.items(),
                        key=lambda item: item[1],
                        reverse=True
                    )
                )
        return json_overview

    def write_closest_canonicals_log_to_file(self):
        json_overview = self.compute_json_overview_dict_for_closest_canonicals()

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
                       self.parser_args.gffcompare_gtf + "\n\n")
            file.write("Reference GTF-file:\n" +
                       self.parser_args.reference_gtf + "\n\n")
            file.write("Reference FASTA-file:\n" +
                       self.parser_args.fasta_path + "\n\n")
            file.write("Specified offset: " +
                       str(self.parser_args.offset) + "\n")
            file.write("Class codes: " +
                       str(self.parser_args.class_codes) + "\n")
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
                           " | " + str(value['closest_canonical']) + " |\n")
            file.close()

        output_manager.output_line({
            "line": "Results written to file: " + FASTA_OVERVIEW_FILE,
            "is_info": True
        })

    def execute_log_file_creation(self):
        output_manager.output_line({
            "line": "Creating log-files",
            "is_info": True
        })
        self.write_closest_canonicals_log_to_file()

        pass
