from services.output_manager import default_output_manager as output_manager


def compute_average_and_sd(findings: dict):
    findings['ins_avg'] = sum(
        [key * value for key, value in findings['insertions'].items()]) / sum(findings['insertions'].values())
    findings['ins_sd'] = (sum(
        [(key - findings['ins_avg']) ** 2 * value for key, value in findings['insertions'].items()]) / sum(
        findings['insertions'].values())) ** 0.5
    findings['del_avg'] = sum(
        [key * value for key, value in findings['deletions'].items()]) / sum(findings['deletions'].values())
    findings['del_sd'] = (sum(
        [(key - findings['del_avg']) ** 2 * value for key, value in findings['deletions'].items()]) / sum(
        findings['deletions'].values())) ** 0.5


def make_prediction(findings: dict, location_type: str):
    total_cases = sum(findings['insertions'].values())
    if total_cases < 5:
        return

    compute_average_and_sd(findings)

    if location_type == "start":
        canonicals = ["AG", "AC"]
    else:
        canonicals = ["GT", "GC", "AT"]

    # if findings['closest_canonical'][1] not in canonicals:
    #     return

    del_max_value = [k for k, v in findings['deletions'].items(
    ) if v == max(findings['deletions'].values())]
    if len(del_max_value) > 1:
        return

    threshold = 0.7
    nucleotides_exceeding_treshold = 0
    for value in findings['del_pos_distr']:
        if value / total_cases > threshold:
            nucleotides_exceeding_treshold += 1

    # if del_max_value[0] == findings['closest_canonical'][2] and findings['closest_canonical'][2] != 0 and nucleotides_exceeding_treshold == del_max_value[0]:
    #     findings['error_detected'] = True

    if del_max_value[0] != 0 and nucleotides_exceeding_treshold == del_max_value[0]:
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


def execute_error_prediction(intron_site_dict: dict):
    output_manager.output_line({
        "line": "PREDICTING ERRORS",
        "is_title": True
    })
    for value in intron_site_dict.values():
        for findings in value["extracted_information"].values():
            make_prediction(findings, value["location_type"])

    count_predicted_errors(intron_site_dict)

    output_manager.output_line({
        "line": "error prediction finished",
        "is_info": True
    })
