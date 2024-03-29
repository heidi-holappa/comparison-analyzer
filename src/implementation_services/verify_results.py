from services.output_manager import default_output_manager as output_manager
from services.log_manager import default_log_manager as log_manager


def verify_results(parser_args, intron_site_dict: dict, matching_cases_dict: dict, class_codes_and_transcripts: dict):
    output_manager.output_line({
        "line": "VERIFYING RESULTS",
        "is_title": True
    })
    results = {
        'TP': {
            'left': {},
            'right': {},
            'closest_canonical_matches': 0
        },
        'FP': {
            'left': {},
            'right': {},
            'closest_canonical_matches': 0
        },
        'unverified_cases': {
            'left': {},
            'right': {},
            'closest_canonical_matches': 0
        },
    }
    verified_cases = 0
    debug_true_positives_dict = {}
    debug_false_positives_dict = {}
    debug_unverified_cases = {}
    unverified_cases_class_codes = {}
    # debug_errors = []
    for key, value in intron_site_dict.items():
        directions = ['right', 'left']
        for direction in directions:
            if value['extracted_information'][direction]['error_detected']:
                deletions = value['extracted_information'][direction]['deletions']
                most_common_del = max(deletions, key=deletions.get)
                matching_case = matching_cases_dict.get(key)
                case_key = key

                closest_canonical_distance = value['extracted_information'][direction]['closest_canonical'][2]
                case_value = {
                    'direction': direction,
                    'strand': value['strand'],
                    'deletions': deletions,
                    'del_pos_distr': value['extracted_information'][direction]['del_pos_distr'],
                    'most_common_del': most_common_del,
                    'most_common_del_pair': value['extracted_information'][direction]['most_common_del_pair'],
                    'del_avg': "{:.2f}".format(value['extracted_information'][direction]['del_avg']),
                    'del_sd': "{:.2f}".format(value['extracted_information'][direction]['del_sd']),
                    'aggressive_strategy': parser_args.no_canonicals,
                }

                if not matching_case:
                    if value['transcript_id'] in class_codes_and_transcripts:
                        class_code = class_codes_and_transcripts[value['transcript_id']]
                        if class_code not in unverified_cases_class_codes:
                            unverified_cases_class_codes[class_code] = 0
                        unverified_cases_class_codes[class_code] += 1
                    debug_unverified_cases[case_key] = case_value
                    # debug_errors.append(str(key) + "\n")
                    if most_common_del not in results['unverified_cases'][direction]:
                        results['unverified_cases'][direction][most_common_del] = 0
                    results['unverified_cases'][direction][most_common_del] += 1
                    if closest_canonical_distance == most_common_del:
                        results['unverified_cases']['closest_canonical_matches'] += 1
                    continue
                verified_cases += 1
                offset = abs(matching_case['offset'])
                # predicted_offset = value['extracted_information'][direction]['closest_canonical'][2]

                if offset == most_common_del:
                    if most_common_del not in results['TP'][direction]:
                        results['TP'][direction][most_common_del] = 0
                    results['TP'][direction][most_common_del] += 1
                    debug_true_positives_dict[case_key] = case_value
                    debug_true_positives_dict[case_key]['offset'] = offset
                    if closest_canonical_distance == offset:
                        results['TP']['closest_canonical_matches'] += 1
                    break
                else:
                    if most_common_del not in results['FP'][direction]:
                        results['FP'][direction][most_common_del] = 0
                    results['FP'][direction][most_common_del] += 1
                    debug_false_positives_dict[case_key] = case_value
                    debug_false_positives_dict[case_key]['offset'] = offset
                    if closest_canonical_distance == offset:
                        results['FP']['closest_canonical_matches'] += 1

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
    output_manager.output_line({
        "line": "Unverified cases: " + str(results['unverified_cases']) + ", total: " + str(sum(results['unverified_cases']['left'].values()) + sum(results['unverified_cases']['right'].values())),
        "is_info": True
    })
    if unverified_cases_class_codes:
        output_manager.output_line({
            "line": "Unverified cases class codes: " + str(unverified_cases_class_codes),
            "is_info": True
        })
    # if debug_errors:
    #     log_manager.debug_logs['verifying_results_keys_not_found'] = debug_errors
    if debug_unverified_cases:
        log_manager.debug_logs['results_unverified_cases'] = debug_unverified_cases
    if debug_true_positives_dict:
        log_manager.debug_logs['results_true_positives'] = debug_true_positives_dict
    if debug_false_positives_dict:
        log_manager.debug_logs['results_false_positives'] = debug_false_positives_dict
