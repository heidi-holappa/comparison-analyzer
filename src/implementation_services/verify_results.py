from services.output_manager import default_output_manager as output_manager
from services.log_manager import default_log_manager as log_manager


def verify_results(intron_site_dict: dict, matching_cases_dict: dict):
    output_manager.output_line({
        "line": "VERIFYING RESULTS",
        "is_title": True
    })
    results = {
        'TP': 0,
        'FP': 0
    }
    debug_errors = []
    for key, value in intron_site_dict.items():
        if value['extracted_information']['right']['error_detected'] or \
                value['extracted_information']['left']['error_detected']:

            case = matching_cases_dict.get(key)
            if not case:
                debug_errors.append(str(key))
                continue
            offset = case['offset']

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
    if debug_errors:
        log_manager.debug_logs['verifying_results_keys_not_found'] = debug_errors
