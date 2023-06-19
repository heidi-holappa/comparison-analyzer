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
        expected_result, expected_debug_list = False, [[]]
        result, debug_list = alignment_parser.count_indels_from_cigar_codes_in_given_window(
            cigar_tuples, aligned_location, loc_type)
        self.assertEqual((result, debug_list),
                         (expected_result, expected_debug_list))

    def test_indels_are_counted_correctly(self):
        cigar_tuples = [(0, 20), (2, 3), (1, 2), (0, 10)]
        self.parser.window_size = 8
        aligned_location = 20
        loc_type = "start"
        expected_errors, expected_debug_list = False, [2, 2, 2, 1, 1, 0, 0, 0]

        errors, debug_list = alignment_parser.count_indels_from_cigar_codes_in_given_window(
            cigar_tuples, aligned_location, loc_type)
        self.assertEqual((errors, debug_list[0]),
                         (expected_errors, expected_debug_list))

    def test_full_window_of_dels_returns_true_for_errors(self):
        cigar_tuples = [(0, 20), (2, 8), (1, 2), (0, 10)]
        self.parser.window_size = 8
        aligned_location = 20
        loc_type = "start"
        expected_errors, expected_debug_list = True, [2, 2, 2, 2, 2, 2, 2, 2]

        errors, debug_list = alignment_parser.count_indels_from_cigar_codes_in_given_window(
            cigar_tuples, aligned_location, loc_type)
        self.assertEqual((errors, debug_list[0]),
                         (expected_errors, expected_debug_list))


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

    def test_read_with_ten_aligned_pairs_loc_type_start_window_size_two_returns_insertion_of_one_for_selected_location(self):
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

    def test_read_with_ten_aligned_pairs_loc_type_end_window_size_two_returns_insertion_of_two_for_selected_location(self):
        """
        List items 3,4 are iterated over and the number of insertions is counted
        """
        aligned_pairs = [(0, None), (1, None), (2, None), (3, None), (4, 105),
                         (5, 105), (6, 106), (7, 107), (8, 108), (9, 109)]
        location = 5
        self.parser.window_size = 2
        loc_type = "end"
        expected_output = {
            "insertions": {
                1: 1
            },
            "deletions": {}
        }
        self.parser.process_read(aligned_pairs, location, loc_type)
        result = self.parser.case_count
        self.assertEqual(result, expected_output)

    def test_read_with_ten_pairs_loc_type_end_with_deletion_counts_correct_amount_of_deletions(self):
        """
        List indexes 2,3,4 are iterated over and the number of deletions is counted        
        """
        aligned_pairs = [(0, 100), (1, 101), (None, 103), (None, 104), (3, 105),
                         (4, 105), (5, 106), (6, 107), (7, 108), (8, 109)]
        location = 4
        self.parser.window_size = 2
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

    def test_read_with_ten_pairs_loc_type_start_with_deletion_counts_correct_amount_of_deletions(self):
        """
        List indexes 2,3,4 are iterated over and the number of deletions is counted        
        """
        aligned_pairs = [(0, 100), (1, 101), (None, 103), (None, 104), (3, 105),
                         (4, 105), (5, 106), (6, 107), (7, 108), (8, 109)]
        location = 2
        self.parser.window_size = 3
        loc_type = "start"
        expected_output = {
            "insertions": {},
            "deletions": {
                2: 1
            }
        }
        self.parser.process_read(aligned_pairs, location, loc_type)
        result = self.parser.case_count
        self.assertEqual(result, expected_output)

    def test_window_size_count_of_deletions_returns_true_as_first_arg(self):
        aligned_pairs = [(0, 100), (1, 101), (None, 103),
                         (None, 104), (3, 105)]
        location = 2
        self.parser.window_size = 2
        loc_type = "start"
        result, debug_list = self.parser.process_read(
            aligned_pairs, location, loc_type)
        expected_output = True
        self.assertEqual(result, expected_output)

    def test_error_log_file_is_written_correctly(self):
        filepath = self.parser.error_file_output_dir
        if os.path.exists(filepath):
            os.remove(filepath)
        errors = ["error1"]
        self.parser.write_alignment_errors_to_file(errors)
        with open(filepath, "r") as f:
            result = f.read()
        print(result)

        self.assertTrue("error1" in result)
        if os.path.exists(filepath):
            os.remove(filepath)


class TestBamReader(TestCase):

    def setUp(self):
        self.parser = AlignmentParser()
        test_bam_file = os.path.join(
            TEST_FILE_DIR, "Mouse.ONT.R9.4.sim.RE.no_gtf.transcript925.ch1.nnic.bam")
        self.parser.initialize_file(test_bam_file)

    def test_bam_reader_runs_without_errors(self):
        reads_and_locations = {
            "ENSMUST00000208994_1011_aligned_5112815_F_38_212_77": [(36796992, "end")]
        }
        self.parser.process_bam_file(reads_and_locations)
