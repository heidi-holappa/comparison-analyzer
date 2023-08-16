import os

from unittest import TestCase
import pickle

from config import LOG_DIR

from services.save_point_handler import get_matching_cases
from services.save_point_handler import get_intron_cases


class TestSavePointHandler(TestCase):

    def setUp(self):
        d = {'a': 1, 'b': 2}
        self.pickle_path = os.path.join(LOG_DIR, 'test.pkl')
        with open(self.pickle_path, 'wb') as f:
            pickle.dump(d, f)

    def test_pickle_file_is_returned_correctly(self):
        d = {'a': 1, 'b': 2}

        imported_data = get_matching_cases(self.pickle_path)

        self.assertEqual(imported_data, d)

    def test_if_pickle_file_is_not_found_none_is_returned_for_matching_cases_dict(self):
        imported_data = get_matching_cases('not_found.pkl')

        self.assertIsNone(imported_data)

    def test_if_pickle_file_is_not_found_none_is_returned_for_intron_cases_dict(self):
        imported_data = get_intron_cases('not_found.pkl')

        self.assertIsNone(imported_data)
