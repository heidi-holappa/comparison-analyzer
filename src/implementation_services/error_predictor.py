
def make_prediction(findings: dict):
    total_cases = sum(findings['insertions'].values())
    if total_cases < 3:
        return
    insertion_average = sum(
        [key * value for key, value in findings['insertions'].items()]) / total_cases
    deletion_average = sum(
        [key * value for key, value in findings['deletions'].items()]) / total_cases
    average_treshold = 2
    if (insertion_average > average_treshold or deletion_average > average_treshold) and findings['closest_canonical'][2] != 0:
        findings['insertion_error'] = True


def execute_error_prediction(intron_site_dict: dict):
    for value in intron_site_dict.values():
        for findings in value["extracted_information"].values():
            make_prediction(findings)
