import os
import json

from services.output_manager import default_output_manager as output_manager
from services.graph_manager import default_graph_manager as graph_manager
from config import LOG_FILE_DIR, FASTA_OVERVIEW_FILE


class LogManager:

    def __init__(self):
        self.matching_cases_dict = {}
        self.debug_logs = {}
        pass

    def compute_indel_results(self):
        indel_results = {}
        for matching_case in self.matching_cases_dict.values():
            for type, count in matching_case['indel_errors'].items():
                key = (
                    type, matching_case['strand'], matching_case['location_type'], matching_case['offset'])
                if key not in indel_results:
                    indel_results[key] = {}
                if count not in indel_results[key]:
                    indel_results[key][count] = 0
                indel_results[key][count] += 1

    def compute_json_overview_dict_for_closest_canonicals(self):
        count = {
            'cases_with_indel_errors': 0,
            'cases_with_indel_errors_not_found': 0,
        }
        json_overview = {}
        for value in self.matching_cases_dict.values():
            strand, location_type = value['strand'], value['location_type']
            offset = value['offset']
            if 'indel_errors' not in value:
                count['cases_with_indel_errors_not_found'] += 1
                continue
            for indel_errors_key in value['indel_errors']:
                count['cases_with_indel_errors'] += 1
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
        print(count)
        return json_overview

    def write_closest_canonicals_log_to_file(self, parser_args):
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
                       parser_args.gffcompare_gtf + "\n\n")
            file.write("Reference GTF-file:\n" +
                       parser_args.reference_gtf + "\n\n")
            file.write("Reference FASTA-file:\n" +
                       parser_args.reference_fasta + "\n\n")
            file.write("Specified offset: " +
                       str(parser_args.offset) + "\n")
            file.write("Class codes: " +
                       str(parser_args.class_codes) + "\n")
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

    def write_debug_files(self):
        for log_name, log_values in self.debug_logs.items():
            filepath = os.path.join(LOG_FILE_DIR, 'debug_' + log_name + '.log')
            with open(filepath, "w") as file:
                for entry_key, entry_values in log_values.items():
                    file.write(f"{entry_key}\t{entry_values}\n")

    def execute_log_file_creation(self, matching_cases_dict: dict, parser_args):

        output_manager.output_line({
            "line": "Creating log-files",
            "is_info": True
        })

        self.matching_cases_dict = matching_cases_dict

        self.write_closest_canonicals_log_to_file(parser_args)
        if parser_args.extended_debug:
            self.write_debug_files()

        output_manager.output_footer()
        output_manager.write_log_file()

        pass


default_log_manager = LogManager()
