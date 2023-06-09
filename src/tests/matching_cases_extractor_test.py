from unittest import TestCase

from tests.sample_file_management import default_test_file_manager as file_manager
from services.offset_computation import execute_offset_computation
from services.extract_matching_cases import MatchingCasesExtractor


class TestMatchingCasesExtractor(TestCase):

    def setUp(self):
        file_manager.initialize_test_files()
        self.reference_db = file_manager.reference_db
        self.gffcompare_db = file_manager.gffcompare_db
        class_code = "j"
        self.offset_results = execute_offset_computation(
            class_code, self.gffcompare_db, self.reference_db)
        self.offset = 4

    def test_extract_matching_cases(self):
        extractor = MatchingCasesExtractor(
            self.offset_results,
            self.offset,
            self.reference_db
        )
        matching_cases_dict = extractor.extract_candidates_matching_selected_offset()
        self.assertEqual(len(matching_cases_dict), 2)

    def tearDown(self) -> None:
        file_manager.remove_test_files()
