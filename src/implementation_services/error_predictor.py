from services.output_manager import default_output_manager as output_manager
from config import PRED_PREC_OF_ALL_CASES_TRESHOLD
from config import PRED_MIN_CASES_THRESHOLD
from config import PRED_ACCEPTED_OFFSET_CASES


def compute_average_and_sd(findings: dict):
    del_sum = sum(findings['deletions'].values())
    ins_sum = sum(findings['insertions'].values())
    findings['ins_avg'] = sum(
        [key * value for key, value in findings['insertions'].items()]) / ins_sum
    findings['ins_sd'] = (sum(
        [(key - findings['ins_avg']) ** 2 * value for key, value in findings['insertions'].items()]) / ins_sum) ** 0.5
    findings['del_avg'] = sum(
        [key * value for key, value in findings['deletions'].items()]) / del_sum
    findings['del_sd'] = (sum(
        [(key - findings['del_avg']) ** 2 * value for key, value in findings['deletions'].items()]) / del_sum) ** 0.5


def verify_sublist_largest_values_exists(lst, n):
    largest_values = set(sorted(lst, reverse=True)[:n])
    count = 0

    for num in lst:
        if num in largest_values:
            count += 1
            if count >= n:
                return True
        else:
            count = 0

    return False


def make_prediction(parser_args, findings: dict, location_type: str, strand: str):

    # Constants
    total_cases_threshold = PRED_MIN_CASES_THRESHOLD
    count_proportion_threshold = PRED_PREC_OF_ALL_CASES_TRESHOLD
    accepted_offset_cases = PRED_ACCEPTED_OFFSET_CASES

    total_cases = sum(findings['insertions'].values())
    suported_strands = ['+', '-']

    if total_cases < total_cases_threshold or strand not in suported_strands:
        return

    compute_average_and_sd(findings)

    # if location_type == "start":
    #     canonicals = ["AG", "AC"]
    # else:
    #     canonicals = ["GT", "GC", "AT"]

    # if findings['closest_canonical'][1] not in canonicals:
    #     return

    del_most_common_case = [k for k, v in findings['deletions'].items(
    ) if v == max(findings['deletions'].values())]

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
    canonicals = possible_canonicals[strand][location_type]

    # If aggressive search is not enabled and the most common deletion pair is not a canonical, do nothing
    if not parser_args.no_canonicals and findings['most_common_del_pair'] not in canonicals:
        return

    # If a distinct most common deletion count does not exists or
    # most common case is not among accepted offset cases, do nothing
    if len(del_most_common_case) > 1 or del_most_common_case[0] not in accepted_offset_cases:
        return

    # Compute count for nucleotides exceeding distribution threshold
    nucleotides_exceeding_treshold = 0
    for value in findings['del_pos_distr']:
        if value / total_cases > count_proportion_threshold:
            nucleotides_exceeding_treshold += 1

    # if del_max_value[0] == findings['closest_canonical'][2] and findings['closest_canonical'][2] != 0 and nucleotides_exceeding_treshold == del_max_value[0]:
    #     findings['error_detected'] = True

    if parser_args.no_canonicals:
        consentration_exists = verify_sublist_largest_values_exists(
            findings['del_pos_distr'], del_most_common_case[0])
        threshold_exceeds = bool(
            nucleotides_exceeding_treshold >= del_most_common_case[0])
        if consentration_exists and threshold_exceeds:
            findings['error_detected'] = True
    elif parser_args.very_conservative:
        consentration_exists = verify_sublist_largest_values_exists(
            findings['del_pos_distr'], del_most_common_case[0])
        threshold_exceeds = bool(
            nucleotides_exceeding_treshold >= del_most_common_case[0])
        canonical_matches = bool(
            findings['closest_canonical'][2] == del_most_common_case[0])
        if consentration_exists and threshold_exceeds and canonical_matches:
            findings['error_detected'] = True
        pass
    else:
        findings['error_detected'] = True


def count_predicted_errors(intron_site_dict: dict):
    count_of_errors = 0
    for value in intron_site_dict.values():
        for findings in value["extracted_information"].values():
            if findings['error_detected']:
                count_of_errors += 1
                break
    output_manager.output_line({
        "line": f"Predicted errors: {count_of_errors}",
        "is_info": True
    })


def execute_error_prediction(parser_args, intron_site_dict: dict):
    output_manager.output_line({
        "line": "PREDICTING ERRORS",
        "is_title": True
    })
    for value in intron_site_dict.values():
        for findings in value["extracted_information"].values():
            make_prediction(parser_args, findings,
                            value["location_type"], value["strand"])

    count_predicted_errors(intron_site_dict)

    output_manager.output_line({
        "line": "error prediction finished",
        "is_info": True
    })
