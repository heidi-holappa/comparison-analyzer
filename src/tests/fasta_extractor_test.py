from unittest import TestCase
import pytest
from services.fasta_extractor import FastaExtractor
from tests.sample_file_management import default_test_file_manager as file_manager


class TestFastaExtractor(TestCase):

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def setUp(self):
        """
            Note: these values have been altered to match the range of the sample data.
            The extracted matches should be 'GT' and 'AG'. This needs to be changed
            when the sample data is updated.
        """
        matching_cases_dict = {
            ('transcript1.chr1.nnic', 210): {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '+',
                'location_type': 'start',
                'exon_number': 3,
                'location': 210,
                'seq_id': 'chr1'
            },
            ('transcript1.chr1.nnic', 45): {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '+',
                'location_type': 'end',
                'exon_number': 4,
                'location': 45,
                'seq_id': 'chr1'
            },
            ('transcript2.chr1.nnic', 150): {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '-',
                'location_type': 'end',
                'exon_number': 3,
                'location': 45,
                'seq_id': 'chr1'
            },
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
        assert 'CLOSEST CANONICALS' in captured.out

    def test_successful_execution_of_extraction_updates_dictionary(self):
        extractor = FastaExtractor(self.fasta_config)
        extractor.execute_fasta_extraction()
        result = bool(
            'closest_canonical' in extractor.matching_cases_dict[('transcript1.chr1.nnic', 210)])
        print(extractor.matching_cases_dict)
        self.assertTrue(result)

    def test_extracting_nucleotides_returns_correct_values(self):
        extractor = FastaExtractor(self.fasta_config)
        extractor.initialize_fasta()
        coordinates = ("chr1", 100, 116)
        nucleotides = extractor.extract_characters_at_given_coordinates(
            coordinates)
        expected_result = 'TTTGTTATCTTCCTGG'
        self.assertEqual(nucleotides, expected_result)

    def test_finding_closest_canonicals_returns_correct_values_when_both_matches_are_found(self):
        extractor = FastaExtractor(self.fasta_config)
        dict_entry = ('transcript1.chr1.nnic', 210)
        nucleotides = 'GAAAGCAAGTATTTTG'
        canonicals = ['GT', 'GC', 'AT']
        extractor.find_closest_canonicals(
            nucleotides, dict_entry, canonicals)
        left = extractor.matching_cases_dict[dict_entry]['closest_canonical']['left']
        right = extractor.matching_cases_dict[dict_entry]['closest_canonical']['right']
        self.assertEqual(left, ("GC", "GT", 4))
        self.assertEqual(right, ("AT", "GT", 2))

    def test_finding_closest_canonicals_returns_correct_values_when_only_right_match_is_found(self):
        extractor = FastaExtractor(self.fasta_config)
        dict_entry = ('transcript1.chr1.nnic', 210)
        nucleotides = 'GAAAACAAGTATTTTG'
        canonicals = ['GT', 'GC', 'AT']
        extractor.find_closest_canonicals(
            nucleotides, dict_entry, canonicals)
        left = extractor.matching_cases_dict[dict_entry]['closest_canonical']['left']
        right = extractor.matching_cases_dict[dict_entry]['closest_canonical']['right']
        self.assertEqual(left, ("GT", "GT", 0))
        self.assertEqual(right, ("AT", "GT", 2))

    def test_finding_closest_canonicals_returns_correct_values_when_only_left_match_is_found(self):
        extractor = FastaExtractor(self.fasta_config)
        dict_entry = ('transcript1.chr1.nnic', 210)
        nucleotides = 'GAAAGCAAGTCTTTTG'
        canonicals = ['GT', 'GC', 'AT']
        extractor.find_closest_canonicals(
            nucleotides, dict_entry, canonicals)
        left = extractor.matching_cases_dict[dict_entry]['closest_canonical']['left']
        right = extractor.matching_cases_dict[dict_entry]['closest_canonical']['right']
        self.assertEqual(left, ("GC", "GT", 4))
        self.assertEqual(right, ("GT", "GT", 0))
