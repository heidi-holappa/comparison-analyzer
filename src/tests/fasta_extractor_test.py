from unittest import TestCase
import pytest
from services.fasta_extractor import FastaExtractor
from tests.sample_file_management import default_test_file_manager as file_manager


class TestFastaExtractor(TestCase):

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def setUp(self):
        matching_cases_dict = {
            ('transcript1.chr1.nnic', 'ENSMUST00000208994.2', '+', 3, 'end'): 36762532,
            ('transcript1.chr1.nnic', 'ENSMUST00000208994.2', '+', 5, 'end'): 36773030
        }
        self.fasta_config = {
            "fasta_path": file_manager.reference_fasta,
            "offset": 4,
            "gffcompare_gtf": file_manager.gffcompare_gtf,
            "reference_gtf": file_manager.reference_gtf,
            "class_codes": "j",
            "matching_cases_dict": matching_cases_dict
        }

    def test_with_no_missing_values_no_errors_are_detected(self):
        extractor = FastaExtractor(self.fasta_config)
        extractor.initialize_fasta()
        errors = extractor.check_errors()
        self.assertEqual(errors, False)

    def test_missing_offset_causes_error(self):
        self.fasta_config['offset'] = -1
        extractor = FastaExtractor(self.fasta_config)
        extractor.initialize_fasta()
        errors = extractor.check_errors()
        self.assertEqual(errors, True)

    def test_empty_mathing_cases_dict_causes_error(self):
        self.fasta_config['matching_cases_dict'] = {}
        extractor = FastaExtractor(self.fasta_config)
        errors = extractor.check_errors()
        self.assertEqual(errors, True)

    def test_an_empty_string_title_generates_a_line_of_selected_character(self):
        extractor = FastaExtractor(self.fasta_config)
        extractor.output_section_header()
        captured = self.capsys.readouterr()
        assert 'FASTA EXTRACTION' in captured.out

    def test_successful_execution_of_extraction_returns_a_dictionary(self):
        extractor = FastaExtractor(self.fasta_config)
        result = extractor.execute_fasta_extraction()
        self.assertIsInstance(result, dict)
