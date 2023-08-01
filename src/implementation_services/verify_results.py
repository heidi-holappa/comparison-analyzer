from services.output_manager import default_output_manager as output_manager
from services.log_manager import default_log_manager as log_manager


def verify_results(intron_site_dict: dict, matching_cases_dict: dict):
    output_manager.output_line({
        "line": "VERIFYING RESULTS",
        "is_title": True
    })
    results = {
        'TP': {
            'left': {},
            'right': {}
        },
        'FP': {
            'left': {},
            'right': {}
        },
    }
    verified_cases = 0
    debug_true_positives_dict = {}
    debug_false_positives_dict = {}
    debug_errors = []
    for key, value in intron_site_dict.items():
        directions = ['right', 'left']
        for direction in directions:
            if value['extracted_information'][direction]['error_detected']:
                deletions = value['extracted_information'][direction]['deletions']
                most_common_del = max(deletions, key=deletions.get)
                case = matching_cases_dict.get(key)
                if not case:
                    debug_errors.append(str(key) + "\n")
                    continue
                verified_cases += 1
                offset = abs(case['offset'])
                # predicted_offset = value['extracted_information'][direction]['closest_canonical'][2]

                if offset == most_common_del:
                    if most_common_del not in results['TP'][direction]:
                        results['TP'][direction][most_common_del] = 0
                    results['TP'][direction][most_common_del] += 1
                    debug_true_positives_dict[key] = value
                    break
                else:

                    if most_common_del not in results['FP'][direction]:
                        results['FP'][direction][most_common_del] = 0
                    results['FP'][direction][most_common_del] += 1
                    debug_false_positives_dict[key] = value

    output_manager.output_line({
        "line": "Verified cases: " + str(verified_cases),
        "is_info": True
    })

    output_manager.output_line({
        "line": "True positives: " + str(results['TP']) + ", total: " + str(sum(results['TP']['left'].values()) + sum(results['TP']['right'].values())),
        "is_info": True
    })
    output_manager.output_line({
        "line": "False positives: " + str(results['FP']) + ", total: " + str(sum(results['FP']['left'].values()) + sum(results['FP']['right'].values())),
        "is_info": True
    })
    if debug_errors:
        log_manager.debug_logs['verifying_results_keys_not_found'] = debug_errors
    if debug_true_positives_dict:
        log_manager.debug_logs['verifying_results_true_positives'] = debug_true_positives_dict
    if debug_false_positives_dict:
        log_manager.debug_logs['verifying_results_false_positives'] = debug_false_positives_dict
