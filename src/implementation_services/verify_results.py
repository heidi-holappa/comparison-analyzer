from services.output_manager import default_output_manager as output_manager
from services.log_manager import default_log_manager as log_manager


def verify_results(intron_site_dict: dict, matching_cases_dict: dict):
    output_manager.output_line({
        "line": "VERIFYING RESULTS",
        "is_title": True
    })
    results = {
        'TP': 0,
        'FP': {},
    }
    debug_errors = []
    for key, value in intron_site_dict.items():
        directions = ['right', 'left']
        for direction in directions:
            if value['extracted_information'][direction]['error_detected']:

                case = matching_cases_dict.get(key)
                if not case:
                    debug_errors.append(str(key))
                    continue
                offset = case['offset']

                if offset != 0:
                    results['TP'] += 1
                else:
                    most_common_del = max(value['extracted_information'][direction]['deletions'],
                                          key=value['extracted_information'][direction]['deletions'].get)
                    if most_common_del not in results['FP']:
                        results['FP'][most_common_del] = 0
                    results['FP'][most_common_del] += 1

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
