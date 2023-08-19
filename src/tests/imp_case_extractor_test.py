from unittest import TestCase
from tests.sample_file_management import default_test_file_manager as file_manager

from implementation_services.case_extractor import CaseExtractor


class ImplementationCaseExtractorTest(TestCase):

    def setUp(self):
        file_manager.initialize_test_files()
        file_manager.initialize_test_files()
        self.isoquant_db = file_manager.isoquant_db
        self.case_extractor = CaseExtractor()

    def test_transcripts_are_extracted_correctly(self):
        set_of_transcripts = self.case_extractor.extract_transcripts(
            self.isoquant_db)
        self.assertEqual(len(set_of_transcripts), 1)

    def test_correct_number_of_cases_are_extracted(self):
        extracted_cases = self.case_extractor.extract_intron_site_locations(
            self.isoquant_db)
        self.assertEqual(len(extracted_cases), 4)

    def test_extracted_cases_have_correct_value_keys(self):
        extracted_cases = self.case_extractor.extract_intron_site_locations(
            self.isoquant_db)
        list_of_value_keys = ['transcript_id', 'strand',
                              'location_type', 'location', 'seq_id', 'extracted_information']

        for case in extracted_cases.values():
            for key in list_of_value_keys:
                self.assertTrue(key in case)
