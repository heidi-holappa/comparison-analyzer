import os

from unittest import TestCase
from services.alignment_parser import default_alignment_parser as alignment_parser
from services.alignment_parser import AlignmentParser

from config import TEST_FILE_DIR


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


class TestIndelCountingFromCigarCodes(TestCase):
    def setUp(self):
        self.parser = AlignmentParser()

    def test_indel_counter_returns_false_and_an_empty_debug_list_for_given_empty_list(self):
        cigar_tuples = []
        aligned_location = 100
        loc_type = "start"
        strand = "+"
        expected_result, expected_debug_list, expected_result = False, [
            []], {'deletions': 0, 'insertions': 0}
        result, debug_list, result = alignment_parser.count_indels_from_cigar_codes_in_given_window(
            cigar_tuples, aligned_location, loc_type)
        self.assertEqual((result, debug_list, result),
                         (expected_result, expected_debug_list, expected_result))

    def test_indels_are_counted_correctly(self):
        cigar_tuples = [(0, 20), (2, 3), (1, 2), (0, 10)]
        self.parser.window_size = 8
        aligned_location = 27
        loc_type = "end"
        strand = "+"
        expected_errors, expected_debug_list, expected_result = False, [
            2, 2, 2, 1, 1, 0, 0, 0], {'deletions': 3, 'insertions': 2}

        errors, debug_list, result = alignment_parser.count_indels_from_cigar_codes_in_given_window(
            cigar_tuples, aligned_location, loc_type)
        self.assertEqual((errors, debug_list[0], result),
                         (expected_errors, expected_debug_list, expected_result))

    def test_full_window_of_dels_returns_true_for_errors(self):
        cigar_tuples = [(0, 20), (2, 8), (1, 2), (0, 10)]
        self.parser.window_size = 8
        aligned_location = 20
        loc_type = "start"
        strand = "+"
        expected_errors, expected_debug_list, expected_result = True, [
            2, 2, 2, 2, 2, 2, 2, 2], {'deletions': 8, 'insertions': 0}

        errors, debug_list, result = alignment_parser.count_indels_from_cigar_codes_in_given_window(
            cigar_tuples, aligned_location, loc_type)
        self.assertEqual((errors, debug_list[0], result),
                         (expected_errors, expected_debug_list, expected_result))
