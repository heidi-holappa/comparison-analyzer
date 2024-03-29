import os
import sys

from collections import defaultdict
from pathlib import Path
import pickle

if len(sys.argv) < 2:
    print(
        'Usage: python3 intron-case-analyzer.py <directory> [file-must-contain-str]]')
    sys.exit(1)

file_contains = ''
if len(sys.argv) == 3:
    file_contains = sys.argv[2]

directory = sys.argv[1]

files = []
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f) and 'intron-cases' in f:
        if file_contains in f:
            files.append(f)


print("Found %d files" % len(files))
max_len = max([len(Path(file).stem) for file in files]) + len(" (24/24)")


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


def normalize_results(cases: dict):
    total_count = sum(cases.values())
    normalized_cases = {}
    for case, count in cases.items():
        normalized = (count / total_count) * 100
        if normalized > 1:
            normalized_cases[case] = round(normalized, 1)
        elif normalized > 0.01:
            normalized_cases[case] = round(normalized, 3)
        else:
            normalized_cases[case] = round(normalized, 6)
    return normalized_cases


results = {}
normalized_results = {}
count = 0
for file in files:
    count += 1
    filename = Path(file).stem.replace('-intron-cases', '')
    output = filename + '.' * (max_len - len(filename)) + \
        '(%d/%d)' % (count, len(files))
    print("Processing file: %s" % output, end='\r')
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

    normalized_results_to_be_included = {}

    for direction, cases in max_del_count.items():
        normalized_cases = normalize_results(cases)
        normalized_results_to_be_included[direction] = normalized_cases

    results[filename] = max_del_count
    normalized_results[filename] = normalized_results_to_be_included

print('\n')

for filename in sorted(normalized_results):
    for direction, cases in normalized_results[filename].items():
        print(filename, direction, {i: dict(cases)[i] for i in sorted(cases)})


for filename, output in results.items():
    for direction, cases in output.items():
        print(filename, direction, {i: dict(cases)[i] for i in sorted(cases)})

with open("table-form.md", "w") as f:
    f.write("| Filename | direction | n | -1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | \n")
    f.write(
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | \n")
    for filename in sorted(results):
        for direction, cases in results[filename].items():
            f.write("| %s | %s | %d | " %
                    (filename, direction, sum(cases.values())))
            for i in range(-1, 9):
                if i not in cases:
                    f.write(" NA |")
                    continue
                f.write(
                    f" {cases[i]} ({normalized_results[filename][direction][i]}%) |")
            f.write("\n")
