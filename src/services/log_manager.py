import os
import json

from services.output_manager import default_output_manager as output_manager
from services.graph_manager import default_graph_manager as graph_manager
from config import LOG_FILE_DIR, FASTA_OVERVIEW_FILE, OFFSET_LOG
from config import PRED_MIN_CASES_THRESHOLD


class LogManager:

    def __init__(self):
        self.matching_cases_dict = {}
        self.debug_logs = {}
        self.alignment_erros = []
        self.alignment_error_log_filepath = os.path.join(
            LOG_FILE_DIR, "debug_alignment_errors.log")
        self.offset_results = {}

    def write_offset_results_to_file(self):
        if not self.offset_results:
            return
        with open(OFFSET_LOG, 'w', encoding="utf-8") as file:
            file.write(
                "transcript_id\treference_id\tclass_code\tstrand\toffsets\n")
            for key, value in self.offset_results.items():
                file.write(
                    f"{key}\t{value['reference_id']}\t{value['class_code']}\t" +
                    f"{value['strand']}\t{value['offsets']}\n")

    def write_alignment_errors_to_file(self):

        with open(self.alignment_error_log_filepath, "w") as file:
            file.write(
                "qname\ttranscripts\tlocation\talign_location\t" +
                "type\tread.reference_start\tread.reference_end\t" +
                "list of alignments\n")
            file.writelines(self.alignment_erros)

    def compute_indel_results(self):
        indel_results = {}
        count_indel_errors = 0
        for matching_case in self.matching_cases_dict.values():
            if 'indel_errors' not in matching_case:
                count_indel_errors += 1
                self.alignment_erros.append(str(matching_case) + "\n")
                continue
            for indel_type, indel_dict in matching_case['indel_errors'].items():
                key = (
                    indel_type, matching_case['strand'],
                    matching_case['location_type'],
                    matching_case['offset'])
                if key not in indel_results:
                    indel_results[key] = {}
                for indel_dict_key, indel_event_count in indel_dict.items():
                    if indel_dict_key not in indel_results[key]:
                        indel_results[key][indel_dict_key] = 0
                    indel_results[key][indel_dict_key] += indel_event_count

        if count_indel_errors:
            output_manager.output_line({
                "line": "Indel results: number of cases without indel errors: " +
                f"{count_indel_errors}",
                "is_error": True
            })
        return indel_results

    def get_canonicals(self, strand, location_type):
        possible_canonicals = {
            '+': {
                'start': ['AG', 'AC'],
                'end': ['GT', 'GC', 'AT']
            },
            '-': {
                'start': ['AC', 'GC', 'AC'],
                'end': ['CT', 'GT']
            }
        }
        return possible_canonicals[strand][location_type]

    def validate_img_should_be_created_for_closest_canonical_dict_entry(self,
                                                                        key: tuple,
                                                                        nucleotide_pair: tuple,
                                                                        cases: dict):
        if sum(cases.values()) < PRED_MIN_CASES_THRESHOLD:
            return False
        # Key[0] is strand, key[1] is location_type
        canonicals = self.get_canonicals(key[0], key[1])
        if nucleotide_pair[1].upper() not in canonicals:
            return False
        return True

    def generate_closest_canonicals_graphs(self, parser_args, dict_of_canonicals: dict):
        output_manager.output_line({
            "line": "Closest canonicals: creating graphs.",
            "is_info": True
        })
        graphs_created = 0

        for key, nucleotides in dict_of_canonicals.items():
            for nucleotide_pair, cases in nucleotides.items():
                if not self.validate_img_should_be_created_for_closest_canonical_dict_entry(
                    key, nucleotide_pair, cases
                ):
                    continue
                case_count = sum(cases.values())

                title = f"{nucleotide_pair}: strand: {key[0]}, exon location: {key[1]}, " + \
                    f"offset: {key[2]}, " + \
                    f"direction of pair: {key[3]}, n of cases: {case_count}"
                filename = "closest-canonicals." + ".strand_" + str(key[0]) + ".exon-loc-" + \
                    str(key[1]) + ".offset-(" + str(key[2]) + ")" + \
                    ".direction-" + \
                    str(key[3]) + "." + str(nucleotide_pair[0]) + \
                    str(nucleotide_pair[1])

                graph_manager.construct_bar_chart_from_dict(
                    graph_values=cases,
                    filename=filename,
                    title=title,
                    x_label=f"Number of errors (n of cases: {case_count})",
                    y_label="Portion of reads",
                )
                graphs_created += 1
        output_manager.output_line({
            "line": f"Closest canonicals: done. A total of {graphs_created} " +
            "graphs created for closest canonicals.",
            "is_info": True
        })

    def write_indel_results_to_file(self, indel_results):
        filepath = os.path.join(LOG_FILE_DIR, 'indel_results.log')
        results = []
        total_reads_in_indel_results = 0
        for key, value in indel_results.items():
            results.append(
                f"in/del: {key[0]}, strand: {key[1]}, exon location: {key[2]}, " +
                f"offset: {key[3]}, n of cases: {sum(value.values())}: {value}\n")
            total_reads_in_indel_results += sum(value.values())

        summary_line = "Indel results: count of reads in indel results: " + \
            f"{total_reads_in_indel_results} " + \
            "(Note: one read can be related to multiple matching cases, " + \
            "or be related to multiple transcripts)."

        with open(filepath, "w", encoding="utf-8") as file:
            file.writelines(results)
            file.write(summary_line)
        output_manager.output_line({
            "line": summary_line,
            "is_info": True
        })

    def validate_indel_grahp_should_be_created(self, parser_args, cases: dict):
        if sum(cases.values()) < int(parser_args.min_reads_for_graph):
            return False
        return True

    def generate_indel_graphs(self, parser_args, indel_results):
        output_manager.output_line({
            "line": "Indel results: generating graphs for insertions and deletions " +
            "found at given locations",
            "is_info": True
        })
        graphs_created = 0
        for key, value in indel_results.items():
            if not self.validate_indel_grahp_should_be_created(parser_args, value):
                continue
            title = f"Type: {key[0]}, strand: {key[1]}, exon location: {key[2]}, " + \
                f"offset: {key[3]}, n of cases: {sum(value.values())}"
            filename = "indel." + str(key[0]) + ".strand_" + str(key[1]) + ".exon-loc-" + \
                str(key[2]) + ".offset-(" + str(key[3]) + ")"

            graph_manager.construct_bar_chart_from_dict(
                graph_values=value,
                filename=filename,
                title=title,
                x_label=f"Number of errors (n of cases: {sum(value.values())})",
                y_label="Portion of reads",
            )
            graphs_created += 1

        output_manager.output_line({
            "line": f"Indel results: graphs for {graphs_created} cases created.",
            "is_info": True
        })

    def compute_closest_canonicals_dict(self):
        closest_canonicals_dict = {}
        for value in self.matching_cases_dict.values():
            strand = value['strand']
            location_type = value['location_type']
            offset = value['offset']

            if 'closest_canonical' not in value:
                continue
            for canonicals_key, canonicals_value in value['closest_canonical'].items():
                dict_key = (strand, location_type, offset, canonicals_key)
                if dict_key not in closest_canonicals_dict:
                    closest_canonicals_dict[dict_key] = {}
                nucleotides, distance = (
                    canonicals_value[0], canonicals_value[1]), canonicals_value[2]
                if nucleotides not in closest_canonicals_dict[dict_key]:
                    closest_canonicals_dict[dict_key][nucleotides] = {}
                if distance not in closest_canonicals_dict[dict_key][nucleotides]:
                    closest_canonicals_dict[dict_key][nucleotides][distance] = 0
                closest_canonicals_dict[dict_key][nucleotides][distance] += 1

        return closest_canonicals_dict

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

            if isinstance(log_values, dict):
                with open(filepath, "w") as file:
                    file.write("{\n")
                    for entry_key, entry_values in log_values.items():
                        file.write(f"{entry_key}: {entry_values},\n")
                    file.write("}\n")
            if isinstance(log_values, list):
                with open(filepath, "w") as file:
                    file.writelines(log_values)
            if isinstance(log_values, set):
                with open(filepath, "w") as file:
                    for entry in log_values:
                        file.write(entry + "\n")
            output_manager.output_line({
                "line": f"{log_name} written to: {filepath}",
                "is_info": True
            })

    def generate_output_for_closest_canonicals(self, parser_args):
        closest_canonicals_dict = self.compute_closest_canonicals_dict()
        self.generate_closest_canonicals_graphs(
            parser_args, closest_canonicals_dict)

    def generate_output_for_indels(self, parser_args):
        indel_results = self.compute_indel_results()
        self.generate_indel_graphs(parser_args, indel_results)
        self.write_indel_results_to_file(indel_results)

    def execute_log_file_creation(self, matching_cases_dict: dict, parser_args):

        output_manager.output_line({
            "line": "CREATING LOG-FILES AND GRAPHS",
            "is_title": True
        })

        self.matching_cases_dict = matching_cases_dict

        self.generate_output_for_closest_canonicals(parser_args)
        self.generate_output_for_indels(parser_args)

        if parser_args.extended_debug:
            self.debug_logs["matching_cases_dict"] = self.matching_cases_dict
            self.write_debug_files()

        pass


default_log_manager = LogManager()
