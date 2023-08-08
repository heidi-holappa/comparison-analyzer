import os
from collections import defaultdict
from pathlib import Path
import pickle

directory = '/home/holaphei/koulutyot/lv-2022-2023/BI-summer-trainee/isoquant-summer-trainee/ca-task-2/comparison-analyzer/saved_results'

files = []
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        if 'intron-cases' in f:
            files.append(f)


def get_intron_cases(intron_save_file: str):
    if os.path.exists(intron_save_file):
        with open(intron_save_file, "rb") as f:
            return pickle.load(f)
    return None


def count_most_common_indel_case(indel_dict: dict):
    most_common_case = [
        k for k, v in indel_dict.items() if v == max(indel_dict.values())]
    if len(most_common_case) > 1 or len(most_common_case) == 0:
        return -1
    return most_common_case[0]


results = {}

for file in files:
    intron_cases = get_intron_cases(file)
    max_del_count = {
        'left': defaultdict(int),
        'right': defaultdict(int)
    }
    if not isinstance(intron_cases, dict):
        continue
    for value in intron_cases.values():
        for direction in ['left', 'right']:
            dels = value['extracted_information'][direction]['deletions']
            most_common_del = count_most_common_indel_case(dels)
            max_del_count[direction][most_common_del] += 1
    filename = Path(file).stem
    results[filename] = max_del_count

for filename, output in results.items():
    for direction, cases in output.items():
        print(filename, direction, dict(cases))
