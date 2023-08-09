from unittest import TestCase
import pytest
from implementation_services.closest_canonicals_extractor import find_closest_canonicals


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
        self.matching_cases_dict = {
            ('transcript1.chr1.nnic', 210): {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '+',
                'location_type': 'start',
                'exon_number': 3,
                'location': 210,
                'seq_id': 'chr1',
                'extracted_information': {
                    'left': {
                        'closest_canonical': ""
                    },
                    'right': {
                        'closest_canonical': ""
                    }
                }
            },
            ('transcript1.chr1.nnic', 45): {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '+',
                'location_type': 'end',
                'exon_number': 4,
                'location': 45,
                'seq_id': 'chr1',
                'extracted_information': {
                    'left': {
                        'closest_canonical': ""
                    },
                    'right': {
                        'closest_canonical': ""
                    }
                }
            },
        }

    def test_finding_closest_canonicals_returns_correct_values_when_both_matches_are_found(self):
        dict_entry = ('transcript1.chr1.nnic', 210)
        nucleotides = 'GAAAGCAAGTATTTTG'
        canonicals = ['GT', 'GC', 'AT']
        find_closest_canonicals(
            nucleotides, dict_entry, canonicals, self.matching_cases_dict)
        left = self.matching_cases_dict[dict_entry]['extracted_information']['left']['closest_canonical']
        right = self.matching_cases_dict[dict_entry]['extracted_information']['right']['closest_canonical']
        self.assertEqual(left, ("GC", "GT", 4))
        self.assertEqual(right, ("AT", "GT", 2))
