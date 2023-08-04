import os
import pickle

from services.output_manager import default_output_manager as output_manager


def get_matching_cases(save_file: str):

    if os.path.exists(save_file):
        output_manager.output_line({
            "line": "Loading matching cases from file: " + save_file,
            "is_info": True
        })
        with open(save_file, "rb") as f:
            return pickle.load(f)
    output_manager.output_line({
        "line": "No matching cases found in file: " + save_file,
        "is_info": True
    })
    return None


def save_matching_cases(save_file: str, matching_cases_dict: dict):
    output_manager.output_line({
        "line": "Saving matching cases to file: " + save_file,
        "is_info": True
    })
    with open(save_file, "wb") as f:
        pickle.dump(matching_cases_dict, f)


def get_intron_cases(intron_save_file: str):
    if os.path.exists(intron_save_file):
        output_manager.output_line({
            "line": "Loading intron cases from file: " + intron_save_file,
            "is_info": True
        })
        with open(intron_save_file, "rb") as f:
            return pickle.load(f)
    output_manager.output_line({
        "line": "No intron cases found in file: " + intron_save_file,
        "is_info": True
    })
    return None


def save_intron_cases(intron_save_file: str, intron_site_dict: dict):
    output_manager.output_line({
        "line": "Saving intron cases to file: " + intron_save_file,
        "is_info": True
    })
    with open(intron_save_file, "wb") as f:
        pickle.dump(intron_site_dict, f)
