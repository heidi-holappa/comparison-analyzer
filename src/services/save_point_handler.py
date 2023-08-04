import os
import pickle


def get_matching_cases(save_file: str):

    if os.path.exists(save_file):
        with open(save_file, "rb") as f:
            return pickle.load(f)
    return None


def save_matching_cases(save_file: str, matching_cases_dict: dict):
    with open(save_file, "wb") as f:
        pickle.dump(matching_cases_dict, f)


def get_intron_cases(intron_save_file: str):
    if os.path.exists(intron_save_file):
        with open(intron_save_file, "rb") as f:
            return pickle.load(f)
    return None


def save_intron_cases(intron_save_file: str, intron_site_dict: dict):
    with open(intron_save_file, "wb") as f:
        pickle.dump(intron_site_dict, f)
