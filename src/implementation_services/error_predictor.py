from services.output_manager import default_output_manager as output_manager


def make_prediction(findings: dict):
    total_cases = sum(findings['insertions'].values())
    if total_cases < 3:
        return

    no_ins_errors = findings['insertions'].get(0, 0)
    no_del_errors = findings['deletions'].get(0, 0)

    # prec_no_ins_errors = no_ins_errors / total_cases
    # prec_no_del_errors = no_del_errors / total_cases
    # error_treshold = 0.25

    if max(findings['deletions'], key=findings['deletions'].get) == findings['closest_canonical'][2] and findings['closest_canonical'][2] != 0:
        findings['error_detected'] = True

    # Let's first focus on deletions
    # if (prec_no_ins_errors < error_treshold or prec_no_del_errors < error_treshold) and findings['closest_canonical'][2] != 0:
    #     findings['error_detected'] = True


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
            make_prediction(findings)

    count_predicted_errors(intron_site_dict)

    output_manager.output_line({
        "line": "error prediction finished",
        "is_info": True
    })
