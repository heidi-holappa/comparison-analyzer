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

    def test_extract_matching_cases_for_offset_four_gives_correct_output(self):
        extractor = MatchingCasesExtractor(
            self.offset_results,
            self.offset,
            self.reference_db
        )
        matching_cases_dict = extractor.extract_candidates_matching_selected_offset()
        self.assertEqual(len(matching_cases_dict), 2)

    def test_extract_matching_cases_for_offset_zero_returns_correct_result(self):
        extractor = MatchingCasesExtractor(
            self.offset_results,
            0,
            self.reference_db
        )
        matching_cases_dict = extractor.extract_candidates_matching_selected_offset()
        self.assertEqual(len(matching_cases_dict), 12)

    def test_offset_results_with_more_diversity_are_handled_correctly(self):
        custom_offset_results = {
            'transcript1.chr1.nnic': {
                'reference_id': 'ENSMUST00000208994.2',
                'strand': '+',
                'class_code': 'j',
                'offsets': [(-1, 0), (0, 0), (0, 4), (0, 0)]
            },
            'transcript2.chr1.nnic': {
                'reference_id': 'ENSMUST00000208994.3',
                'strand': '+',
                'class_code': 'j',
                'offsets': [(-1, 0), (0, 0), (4, 0), (0, 0)]
            },
            'transcript3.chr1.nnic': {
                'reference_id': 'ENSMUST00000208994.4',
                'strand': '-',
                'class_code': 'j',
                'offsets': [(-1, 0), (0, 0), (0, 4), (0, 0)]
            },
            'transcript4.chr1.nnic': {
                'reference_id': 'ENSMUST00000208994.5',
                'strand': '-',
                'class_code': 'j',
                'offsets': [(-1, 0), (0, 0), (4, 0), (0, 0)]
            },
            'transcript5.chr1.nnic': {
                'reference_id': 'ENSMUST00000208994.6',
                'strand': '+',
                'class_code': 'j',
                'offsets': [(-1, 0), (0, 0), (0, -4), (0, 0)]
            },
            'transcript6.chr1.nnic': {
                'reference_id': 'ENSMUST00000208994.7',
                'strand': '+',
                'class_code': 'j',
                'offsets': [(-1, 0), (0, 0), (-4, 0), (0, 0)]
            },
            'transcript7.chr1.nnic': {
                'reference_id': 'ENSMUST00000208994.8',
                'strand': '-',
                'class_code': 'j',
                'offsets': [(-1, 0), (0, 0), (0, -4), (0, 0)]
            },
            'transcript8.chr1.nnic': {
                'reference_id': 'ENSMUST00000208994.9',
                'strand': '-',
                'class_code': 'j',
                'offsets': [(-1, 0), (0, 0), (-4, 0), (0, 0)]
            }
        }

        extractor = MatchingCasesExtractor(
            custom_offset_results,
            4,
            self.reference_db
        )

        matching_cases_dict = extractor.extract_candidates_matching_selected_offset()
        print(matching_cases_dict)
        self.assertEqual(len(matching_cases_dict), 1)

    def tearDown(self) -> None:
        file_manager.remove_test_files()
