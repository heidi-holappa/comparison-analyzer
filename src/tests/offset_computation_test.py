import os
from unittest import TestCase
from services.offset_computation import compute_offsets
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


class TestOffsetComputationExecution(TestCase):

    def setUp(self):
        self.test_file_manager = default_test_file_manager
        self.test_file_manager.initialize_test_files()

    def test_execute_offset_computation_returns_a_dictionary_of_results(self):
        class_code = "j"
        results = execute_offset_computation(
            class_code,
            self.test_file_manager.gffcompare_db,
            self.test_file_manager.reference_db,
            False
        )
        self.assertEqual(len(results), 1)

    def tearDown(self):
        self.test_file_manager.remove_test_files()
