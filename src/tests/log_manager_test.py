import os
from unittest import TestCase

from services.log_manager import default_log_manager as log_manager


class TestLogManager(TestCase):

    def setUp(self):
        pass

    def test_indel_results_are_correctly_computated(self):
        log_manager.matching_cases_dict = {
            "case1": {
                "strand": "+",
                "offset": 1,
                "location_type": "start",
                "indel_errors": {
                    "insertions": 2,
                    "deletions": 1
                }
            },
            "case2": {
                "strand": "+",
                "offset": 1,
                "location_type": "start",
                "indel_errors": {
                    "insertions": 3,
                    "deletions": 1
                }
            },
        }
        expected_result = {
            ("insertions", "+", "start", 1): {2: 1, 3: 1},
            ("deletions", "+", "start", 1): {1: 2},
        }

        self.assertEqual(log_manager.compute_indel_results(), expected_result)

    def test_closest_canonicals_are_correctly_computated(self):
        log_manager.matching_cases_dict = {
            "case1": {
                "strand": "+",
                "offset": 1,
                "location_type": "start",
                "indel_errors": {
                    "insertions": 2,
                    "deletions": 1
                },
                "closest_canonical": {
                    "left": ("AT", "GT"),
                    "right": ("GC", "GT")
                }
            },
            "case2": {
                "strand": "+",
                "offset": 1,
                "location_type": "start",
                "indel_errors": {
                    "insertions": 3,
                    "deletions": 1
                },
                "closest_canonical": {
                    "left": ("AT", "GT"),
                    "right": ("GT", "GT")
                }
            },
        }

        expected_result = {
            ('+', 'start', 1, 'left'): {
                ("AT", "GT"): 2
            },
            ('+', 'start', 1, 'right'): {
                ("GC", "GT"): 1,
                ("GT", "GT"): 1
            }
        }
        self.assertEqual(
            log_manager.compute_closest_canonicals_dict(), expected_result)

    def test_errors_are_written_to_file(self):
        if os.path.exists(log_manager.error_file_output_dir):
            os.remove(log_manager.error_file_output_dir)
        log_manager.alignment_erros.append("test_error")
        log_manager.write_alignment_errors_to_file()
        self.assertTrue(os.path.exists(log_manager.error_file_output_dir))

    def tearDown(self):
        if os.path.exists(log_manager.error_file_output_dir):
            os.remove(log_manager.error_file_output_dir)
