import os
from pathlib import Path
import pickle

from config import LOG_DIR


def get_matching_cases(save_file: str):

    if os.path.exists(save_file):
        with open(save_file, "rb") as f:
            return pickle.load(f)
    return None


def save_matching_cases(save_file: str, matching_cases_dict: dict):
    with open(save_file, "wb") as f:
        pickle.dump(matching_cases_dict, f)
