import os
import json

from services.output_manager import default_output_manager as output_manager
from services.graph_manager import default_graph_manager as graph_manager
from config import LOG_FILE_DIR, FASTA_OVERVIEW_FILE, OFFSET_LOG


class LogManager:

    def __init__(self):
        self.matching_cases_dict = {}
        self.debug_logs = {}
        self.alignment_erros = []
        self.alignment_error_log_filepath = os.path.join(
            LOG_FILE_DIR, "debug_alignment_errors.log")
        self.offset_results = {}

    def write_offset_results_to_file(self):
        with open(OFFSET_LOG, 'w', encoding="utf-8") as file:
            file.write(
                "transcript_id\treference_id\tclass_code\tstrand\toffsets\n")
            for key, value in self.offset_results.items():
                file.write(
                    f"{key}\t{value['reference_id']}\t{value['class_code']}\t{value['strand']}\t{value['offsets']}\n")

    def write_alignment_errors_to_file(self):

        with open(self.alignment_error_log_filepath, "w") as file:
            file.write(
                "qname\ttranscripts\tlocation\talign_location\t" +
                "type\tread.reference_start\tread.reference_end\t" +
                "list of alignments\n")
            file.writelines(self.alignment_erros)

    def compute_indel_results(self):
        indel_results = {}
        count_no_indel_errors = 0
        for matching_case in self.matching_cases_dict.values():
            if 'indel_errors' not in matching_case:
                count_no_indel_errors += 1
                self.alignment_erros.append(str(matching_case) + "\n")
                continue
            for type, count in matching_case['indel_errors'].items():
                key = (
                    type, matching_case['strand'], matching_case['location_type'], matching_case['offset'])
                if key not in indel_results:
                    indel_results[key] = {}
                if count not in indel_results[key]:
                    indel_results[key][count] = 0
                indel_results[key][count] += 1
        output_manager.output_line({
            "line": f"Number of cases without indel errors: {count_no_indel_errors}",
            "is_error": True
        })
        return indel_results

    def generate_graphs(self):
        indel_results = self.compute_indel_results()
        output_manager.output_line({
            "line": "Insertions and deletions found at given locations",
            "is_info": True
        })
        total_reads_in_indel_results = 0
        for key, value in indel_results.items():
            title = f"Type: {key[0]}, strand: {key[1]}, exon location: {key[2]}, offset: {key[3]}, n of cases: {sum(value.values())}"
            filename = str(key[0]) + ".strand_" + str(key[1]) + ".exon-loc-" + \
                str(key[2]) + ".offset-(" + str(key[3]) + ")"
            output_manager.output_line({
                "line": f"in/del: {key[0]}, strand: {key[1]}, exon location: {key[2]}, offset: {key[3]}, n of cases: {sum(value.values())}: {value}",
                "is_info": True
            })
            total_reads_in_indel_results += sum(value.values())
            graph_manager.construct_bar_chart_from_dict(
                graph_values=value,
                filename=filename,
                title=title,
                x_label=f"Number of errors (n of cases: {sum(value.values())})",
                y_label="Portion of reads",
            )
        output_manager.output_line({
            "line": f"Total number of reads in indel results: {total_reads_in_indel_results}",
            "is_info": True
        })
        output_manager.output_line({
            "line": f"Graphs for {len(indel_results)} cases created.",
            "is_info": True
        })

    def compute_closest_canonicals_dict(self):
        closest_canonicals_dict = {}
        for value in self.matching_cases_dict.values():
            strand, location_type = value['strand'], value['location_type']
            offset = value['offset']
            for canonicals_key, canonicals_value in value['closest_canonical'].items():
                dict_key = (strand, location_type, offset, canonicals_key)
                canonicals_value_as_str = canonicals_value
                if dict_key not in closest_canonicals_dict:
                    closest_canonicals_dict[dict_key] = {}
                if canonicals_value_as_str not in closest_canonicals_dict[dict_key]:
                    closest_canonicals_dict[dict_key][canonicals_value_as_str] = 0
                closest_canonicals_dict[dict_key][canonicals_value_as_str] += 1

        # for site_locations in closest_canonicals_dict.values():
        #     for results in site_locations.values():
        #         results = dict(
        #             sorted(
        #                 results.items(),
        #                 key=lambda item: item[1],
        #                 reverse=True
        #             )
        #         )

        return closest_canonicals_dict

    def generate_json_overview_dict_for_closest_canonicals(self):
        closest_canonicals_dict = self.compute_closest_canonicals_dict()
        json_overview = {}
        for key, canonical_values in closest_canonicals_dict.items():
            json_overview[str(key)] = {}
            for item in canonical_values:
                json_overview[str(key)][str(item)] = str(
                    canonical_values[item])
        return json_overview

    def write_closest_canonicals_log_to_file(self, parser_args):
        json_overview = self.generate_json_overview_dict_for_closest_canonicals()

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
                       str(parser_args.gffcompare_gtf) + "\n\n")
            file.write("Reference GTF-file:\n" +
                       str(parser_args.reference_gtf) + "\n\n")
            file.write("Reference FASTA-file:\n" +
                       str(parser_args.reference_fasta) + "\n\n")
            file.write("Specified offset: " +
                       str(parser_args.offset) + "\n")
            file.write("Class codes: " +
                       str(parser_args.class_code) + "\n")
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
        output_manager.output_line({
            "line": "CREATING DEBUG LOGS",
            "is_title": True
        })
        self.write_offset_results_to_file()

        output_manager.output_line({
            "line": f"offset results written to: {OFFSET_LOG}",
            "is_info": True
        })

        if self.alignment_erros:
            self.write_alignment_errors_to_file()

        output_manager.output_line({
            "line": f"alignment errors written to: {self.alignment_error_log_filepath}",
            "is_info": True
        })

        for log_name, log_values in self.debug_logs.items():
            filepath = os.path.join(LOG_FILE_DIR, 'debug_' + log_name + '.log')
            with open(filepath, "w") as file:
                for entry_key, entry_values in log_values.items():
                    file.write(f"{entry_key}\t{entry_values}\n")
            output_manager.output_line({
                "line": f"{log_name} written to: {filepath}",
                "is_info": True
            })

    def execute_log_file_creation(self, matching_cases_dict: dict, parser_args):

        output_manager.output_line({
            "line": "Creating log-files",
            "is_info": True
        })

        self.matching_cases_dict = matching_cases_dict

        self.write_closest_canonicals_log_to_file(parser_args)
        self.generate_graphs()

        if parser_args.extended_debug:
            self.write_debug_files()
            self.debug_logs["matching_cases_dict"] = self.matching_cases_dict

        output_manager.output_footer()
        output_manager.write_log_file()

        pass


default_log_manager = LogManager()
