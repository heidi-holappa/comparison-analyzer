import os
import sys

from collections import defaultdict
from pathlib import Path
import pickle

if len(sys.argv) < 2:
    print('Please provide a directory to search for intron cases')
    sys.exit(1)

directory = sys.argv[1]

files = []
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        if 'intron-cases' in f:
            files.append(f)

print("Found %d files" % len(files))
max_len = max([len(Path(file).stem) for file in files])


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
    filename = Path(file).stem
    filename = filename + ' ' * (max_len - len(filename))
    print("Processing file: %s" % filename, end='\r')
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

    results[filename] = max_del_count

print('\n')

for filename, output in results.items():
    for direction, cases in output.items():
        print(filename, direction, sorted(
            dict(cases).items(), key=lambda x: x[1], reverse=True))
