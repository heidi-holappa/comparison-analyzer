from unittest import TestCase
from services.alignment_parser import default_alignment_parser as alignment_parser
from services.alignment_parser import AlignmentParser


class TestCigarParser(TestCase):

    def setUp(self):
        pass

    def test_cigar_string_with_soft_clip_and_one_match_is_parsed_correctly(self):
        cigar = [(4, 50), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 105
        expected_output = 55
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_cigar_string_with_soft_clip_insertion_and_one_match_is_parsed_correctly(self):
        cigar = [(4, 50), (1, 10), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 105
        expected_output = 65
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_cigar_str_with_s_d_i_m_gives_correct_output(self):
        cigar = [(4, 50), (2, 10), (1, 10), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 115
        expected_output = 75
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_cigar_str_with_s_d_n_m_gives_correct_output(self):
        cigar = [(4, 50), (2, 10), (3, 100), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 215
        expected_output = 165
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_cigar_str_with_s_m_i_n_m_gives_correct_output(self):
        cigar = [(4, 50), (0, 10), (1, 10), (3, 100), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 215
        expected_output = 175
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_location_outside_of_cigar_str_returns_minus_one(self):
        cigar = [(4, 50), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 199
        expected_output = -1
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_more_complicated_test_returns_correct_position(self):
        cigar_tuples = [(4, 156), (0, 12), (2, 3), (0, 2), (2, 2), (0, 10), (2, 2), (0, 4), (2, 3), (0, 7), (1, 1), (0, 16), (1, 4), (0, 23), (1, 1), (0, 7),
                        (1, 1), (0, 9), (2, 1), (0, 13), (2, 1), (0, 15), (2, 2), (0, 3), (1, 2), (0, 19), (2, 2), (0, 20), (2, 1), (0, 32), (3, 294), (0, 36), (4, 25)]
        reference_start = 72822568
        reference_end = 73822568
        position = 72823071
        expected_output = 668
        result = alignment_parser.extract_location_from_cigar_string(
            cigar_tuples, reference_start, reference_end, position)
        self.assertEqual(result, expected_output)

    def test_case_that_does_not_consume_any_reference_returns_the_correct_location(self):
        cigar = [(4, 50), (0, 10)]
        reference_start = 100
        reference_end = 160
        location = 100
        expected_output = 50
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_case_that_has_no_reference_consuming_codes_returns_minus_one_as_error(self):
        cigar = [(4, 50), (1, 10)]
        reference_start = 100
        reference_end = 160
        location = 100
        expected_output = -1
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_case_that_has_no_reference_consuming_codes_at_the_end_returns_minus_one_as_error(self):
        cigar = [(4, 50), (0, 10), (1, 10)]
        reference_start = 100
        reference_end = 160
        location = 110
        expected_output = -1
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)

    def test_case_that_has_it_s_location_at_final_match_returns_correct_value(self):
        cigar = [(4, 50), (0, 10), (1, 10)]
        reference_start = 100
        reference_end = 110
        location = 110
        expected_output = 60
        result = alignment_parser.extract_location_from_cigar_string(
            cigar, reference_start, reference_end, location)
        self.assertEqual(result, expected_output)


class TestReadParser(TestCase):

    def setUp(self):
        self.parser = AlignmentParser()

    def test_read_with_no_aligned_pairs_returns_false(self):
        read = []
        location = 100
        loc_type = "start"
        expected_output = False
        result, list = alignment_parser.process_read(read, location, loc_type)
        self.assertEqual(result, expected_output)

    def test_read_with_no_aligned_pairs_does_not_change_the_dictionary_of_counts(self):
        dictionary_of_counts = alignment_parser.case_count
        read = []
        location = 100
        loc_type = "start"
        alignment_parser.process_read(read, location, loc_type)
        self.assertEqual(dictionary_of_counts, alignment_parser.case_count)

    def test_read_with_ten_aligned_pairs_window_size_two_returns_insertion_of_one_for_selected_location(self):
        """
        List items 3,4 are iterated over and the number of insertions is counted
        """
        aligned_pairs = [(0, None), (1, None), (2, None), (3, None), (4, 105),
                         (5, 105), (6, 106), (7, 107), (8, 108), (9, 109)]
        location = 3
        self.parser.window_size = 2
        loc_type = "start"
        expected_output = {
            "insertions": {
                1: 1
            },
            "deletions": {}
        }
        self.parser.process_read(aligned_pairs, location, loc_type)
        result = self.parser.case_count
        self.assertEqual(result, expected_output)

    def test_read_with_ten_pairs_with_deletion_counts_correct_amount_of_deletions(self):
        """
        List indexes 2,3,4 are iterated over and the number of deletions is counted        
        """
        aligned_pairs = [(0, 100), (1, 101), (None, 103), (None, 104), (3, 105),
                         (4, 105), (5, 106), (6, 107), (7, 108), (8, 109)]
        location = 5
        self.parser.window_size = 3
        loc_type = "end"
        expected_output = {
            "insertions": {},
            "deletions": {
                2: 1
            }
        }
        self.parser.process_read(aligned_pairs, location, loc_type)
        result = self.parser.case_count
        self.assertEqual(result, expected_output)
