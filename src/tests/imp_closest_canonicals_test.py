from unittest import TestCase
import pytest

from .sample_file_management import SampleFileManagement

from implementation_services.closest_canonicals_extractor import find_closest_canonicals
from implementation_services.closest_canonicals_extractor import extract_characters_at_given_coordinates
from implementation_services.closest_canonicals_extractor import initialize_fasta
from implementation_services.closest_canonicals_extractor import count_most_common_indel_case
from implementation_services.closest_canonicals_extractor import iterate_intron_sites
from implementation_services.closest_canonicals_extractor import extract_nucleotides_from_most_common_del_location


class TestClosestCanonicalsExtractor(TestCase):

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def setUp(self):
        """
            Note: these values have been altered to match the range of the sample data.
            The extracted matches should be 'GT' and 'AG'. This needs to be changed
            when the sample data is updated.
        """
        self.fasta_path = SampleFileManagement().reference_fasta
        self.intron_site_dict = {
            ('transcript1.chr1.nnic', 210): {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '+',
                'location_type': 'start',
                'exon_number': 3,
                'location': 210,
                'seq_id': 'chr1',
                'extracted_information': {
                    'left': {
                        'closest_canonical': "",
                        'deletions': {0: 1, 1: 2, 4: 10}
                    },
                    'right': {
                        'closest_canonical': "",
                        'deletions': {0: 1, 1: 2, 2: 10}
                    }
                }
            },
            ('transcript1.chr1.nnic', 45): {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': 'X',
                'location_type': 'end',
                'exon_number': 4,
                'location': 45,
                'seq_id': 'chr1',
                'extracted_information': {
                    'left': {
                        'closest_canonical': "",
                        'deletions': {0: 1, 1: 2, 4: 10}
                    },
                    'right': {
                        'closest_canonical': "",
                        'deletions': {0: 1, 1: 2, 4: 10}
                    }
                }
            },
        }

    def test_finding_closest_canonicals_returns_correct_values_when_both_matches_are_found(self):
        dict_entry = ('transcript1.chr1.nnic', 210)
        nucleotides = 'GAAAGCAAGTATTTTG'
        canonicals = ['GT', 'GC', 'AT']
        find_closest_canonicals(
            nucleotides, dict_entry, canonicals, self.intron_site_dict)
        left = self.intron_site_dict[dict_entry]['extracted_information']['left']['closest_canonical']
        right = self.intron_site_dict[dict_entry]['extracted_information']['right']['closest_canonical']
        self.assertEqual(left, ("GC", "GT", 4))
        self.assertEqual(right, ("AT", "GT", 2))

    def test_fasta_extraction_returns_xx_in_case_of_key_error(self):
        fasta = initialize_fasta(self.fasta_path)

        result = extract_characters_at_given_coordinates(
            ('chr3', 10, 12), -1, fasta)
        self.assertEqual(result, 'XX')

    def test_distinct_most_common_case_is_returned(self):
        cases = {0: 10, 1: 2, 3: 0, 4: 20, 5: 1}
        result = count_most_common_indel_case(cases)
        self.assertEqual(result, 4)

    def test_cases_where_strand_is_neither_pos_or_neg_are_omitted(self):
        dict_entry = ('transcript1.chr1.nnic', 45)
        fasta = initialize_fasta(self.fasta_path)
        iterate_intron_sites(self.intron_site_dict, 8, -1, fasta)

        left = self.intron_site_dict[dict_entry]['extracted_information']['left']['closest_canonical']
        right = self.intron_site_dict[dict_entry]['extracted_information']['right']['closest_canonical']
        print(self.intron_site_dict)
        self.assertEqual(left, ("XX", "XX", 0))
        self.assertEqual(right, ("XX", "XX", 0))

    def test_correct_nucleotides_are_extracted(self):
        dict_entry = ('transcript1.chr1.nnic', 210)
        fasta = initialize_fasta(self.fasta_path)
        extract_nucleotides_from_most_common_del_location(
            self.intron_site_dict[dict_entry], fasta)
        left = self.intron_site_dict[dict_entry]['extracted_information']['left']['most_common_case_nucleotides']
        right = self.intron_site_dict[dict_entry]['extracted_information']['right']['most_common_case_nucleotides']
        self.assertEqual(left, 'AG')
        self.assertEqual(right, 'AT')
