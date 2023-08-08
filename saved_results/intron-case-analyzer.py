import os
import pickle

directory = 'saved_results'

files = []
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        if 'intron_cases' in f:
            files.append(f)


def get_intron_cases(intron_save_file: str):
    if os.path.exists(intron_save_file):
        with open(intron_save_file, "rb") as f:
            return pickle.load(f)
    return None


for file in files:
    intron_cases = get_intron_cases(file)
    max_del_count = {
        'left': {},
        'right': {}
    }
    for value in intron_cases.values():
        for 
