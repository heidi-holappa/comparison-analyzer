import os
from unittest import TestCase
from services.offset_computation import compute_offsets
from services.offset_computation import initialize_output_file, write_to_output_file
from services.offset_computation import execute_offset_computation
from config import OFFSET_LOG, TEST_FILE_DIR
from tests.sample_file_management import default_test_file_manager


class TestOffsetComputation(TestCase):

    def setUp(self):
        pass

    def test_offset_for_two_lists_with_one_element_each_should_return_offset_for_those_elements(self):
        ref = [(10, 20)]
        aligned = [(15, 23)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(result, [(5, 3)])

    def test_offset_for_two_references_and_one_aligned_exon_which_is_mapped_to_the_first_reference_exon_is_computed_correctly(self):
        ref = [(10, 20), (30, 40)]
        aligned = [(15, 23)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(result, [(5, 3), (float('-inf'), float('-inf'))])

    def test_offset_for_two_references_and_one_aligned_exon_which_is_mapped_to_the_second_reference_exon_is_computed_correctly(self):
        ref = [(10, 20), (30, 40)]
        aligned = [(25, 33)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(result, [(float('-inf'), float('-inf')), (-5, -7)])

    def test_three_ref_one_aligned_which_is_mapped_to_second_ref_is_computed_correctly(self):
        ref = [(10, 20), (30, 40), (50, 60)]
        aligned = [(25, 33)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(result, [
                         (float('-inf'), float('-inf')), (-5, -7), (float('-inf'), float('-inf'))])

    def test_five_ref_one_aligned_mapped_to_third_ref_is_computed_correctly(self):
        ref = [(10, 20), (30, 40), (50, 60), (70, 80), (90, 100)]
        aligned = [(55, 63)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(result, [(float('-inf'), float('-inf')), (float('-inf'), float(
            '-inf')), (5, 3), (float('-inf'), float('-inf')), (float('-inf'), float('-inf'))])

    def test_offset_for_one_reference_and_two_aligned_exons_for_which_reference_is_mapped_to_first_exon_is_computed_correctly(self):
        ref = [(10, 20)]
        aligned = [(15, 23), (25, 33)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(result, [(5, 3), (float('inf'), float('inf'))])

    def test_offset_for_one_reference_and_five_aligned_exons_for_which_ref_is_mapped_to_third_exon_is_computed_correctly(self):
        ref = [(30, 38)]
        aligned = [(15, 23), (25, 33), (35, 40), (45, 53), (55, 63)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(result, [(float('inf'), float('inf')), (float('inf'), float(
            'inf')), (5, 2), (float('inf'), float('inf')), (float('inf'), float('inf'))])

    def test_five_refs_and_four_offsets_are_mapped_correctly(self):
        ref = [(10, 20), (30, 40), (50, 60), (70, 80), (90, 100)]
        aligned = [(15, 23), (25, 33), (35, 40), (45, 53)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(result, [(5, 3), (float('inf'), float(
            'inf')), (5, 0), (-5, -7), (float('-inf'), float('-inf')), (float('-inf'), float('-inf'))])

    def test_with_a_longer_chain_optimal_matches_are_made(self):
        ref = [(40, 50), (60, 62)]
        aligned = [(10, 20), (50, 60), (61, 63)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(
            result, [(float('inf'), float('inf')), (10, 10), (1, 1)])

    def test_if_no_more_aligned_exons_rest_of_reference_is_mapped_to_minus_infinity(self):
        ref = [(10, 20), (30, 40)]
        aligned = [(15, 23)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(result, [(5, 3), (float('-inf'), float('-inf'))])

    def test_an_alignment_exon_with_no_pair_is_mapped_to_minus_infinity(self):
        aligned = [(10, 20), (30, 40), (50, 60)]
        ref = [(15, 23), (25, 33), (31, 41)]
        result = compute_offsets(aligned, ref)
        self.assertEqual(
            result, [(-5, -3), (float('-inf'), float('-inf')), (-1, -1), (float('inf'), float('inf'))])


class TestOffsetComputationFileManagement(TestCase):

    def test_a_test_file_is_initialized(self):
        self.tearDown()
        initialize_output_file()
        self.assertTrue(os.path.exists(OFFSET_LOG))

    def test_initialized_log_file_only_has_header(self):
        self.tearDown()
        initialize_output_file()
        with open(OFFSET_LOG, 'r') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 1)

    def test_lines_are_written_to_log_file(self):
        self.tearDown()
        initialize_output_file()
        class_code = "j"
        class_code_results_dict = {
            "transcript_1": [(1, 2), (3, 4)],
            "transcript_2": [(5, 6), (7, 8)],
            "transcript_3": [(9, 10), (11, 12)]
        }

        write_to_output_file(class_code_results_dict, class_code)
        with open(OFFSET_LOG, 'r') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 6)

    def tearDown(self) -> None:
        if os.path.exists(OFFSET_LOG):
            os.remove(OFFSET_LOG)


class TestOffsetComputationExecution(TestCase):

    def setUp(self):
        self.test_file_manager = default_test_file_manager
        self.test_file_manager.initialize_test_files()

    def test_execute_offset_computation_returns_a_dictionary_of_results(self):
        class_code = "j"
        results = execute_offset_computation(
            class_code,
            self.test_file_manager.gffcompare_db,
            self.test_file_manager.reference_db
        )
        self.assertEqual(len(results), 1)

    def tearDown(self):
        self.test_file_manager.remove_test_files()
