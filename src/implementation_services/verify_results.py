from services.output_manager import default_output_manager as output_manager


def verify_results(intron_site_dict: dict, matching_cases_dict: dict):
    output_manager.output_line({
        "line": "VERIFYING RESULTS",
        "is_title": True
    })
    results = {
        'TP': 0,
        'FP': 0
    }
    for key, values in intron_site_dict.values():
        if values['extracted_information']['right']['error_found'] == True or \
                values['extracted_information']['left']['error_found'] == True:
            offset = matching_cases_dict[key]['offset']
            if offset != 0:
                results['TP'] += 1
            else:
                results['FP'] += 1

    output_manager.output_line({
        "line": "True positives: " + str(results['TP']),
        "is_info": True
    })
    output_manager.output_line({
        "line": "False positives: " + str(results['FP']),
        "is_info": True
    })
