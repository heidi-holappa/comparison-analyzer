import os
import sys

import pytest
from unittest import TestCase

from services.argument_parser import init_argparser
from services.log_manager import LogManager
from services.log_manager import default_log_manager as log_manager
from config import OFFSET_LOG
from config import PRED_MIN_CASES_THRESHOLD
from config import PRED_PREC_OF_ALL_CASES_TRESHOLD
from config import PRED_ACCEPTED_OFFSET_CASES


class TestLogManager(TestCase):

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def setUp(self):
        sys.argv = [
            'compAna.py',
            '-g=gffcompare.gtf',
            '-r=reference.gtf',
            '-b=file.bam',
            '-a=file.fasta',
            '-t=file.tsv',
            '-o=0 10',
            '-f',
            '-s',
            '-c=j k',
            '-n=true',
            '-v=true',
        ]
        self.parser = init_argparser()

    def test_indel_results_are_correctly_computated(self):
        log_manager.matching_cases_dict = {
            "case1": {
                "strand": "+",
                "offset": 1,
                "location_type": "start",
                "indel_errors": {
                    "insertions": {2: 1, 3: 1},
                    "deletions": {1: 2}
                }
            },
            "case2": {
                "strand": "+",
                "offset": 1,
                "location_type": "start",
                "indel_errors": {
                    "insertions": {3: 1, 4: 2},
                    "deletions": {1: 2, 2: 1}
                }
            },
        }
        expected_result = {
            ("insertions", "+", "start", 1): {2: 1, 3: 2, 4: 2},
            ("deletions", "+", "start", 1): {1: 4, 2: 1},
        }

        self.assertEqual(log_manager.compute_indel_results(), expected_result)

    def test_dictionary_entries_with_no_indel_errrs_are_correctly_computated(self):
        log_manager.matching_cases_dict = {
            "case1": {
                "strand": "+",
                "offset": 1,
                "location_type": "start",
                "closest_canonical": {
                    "left": ("AT", "GT", 2),
                    "right": ("GC", "GT", 4)
                }
            },
            "case2": {
                "strand": "+",
                "offset": 1,
                "location_type": "start",
                "closest_canonical": {
                    "left": ("AT", "GT", 2),
                    "right": ("GT", "GT", 0)
                }
            },
        }
        log_manager.compute_indel_results()
        captured = self.capsys.readouterr()
        assert "2" in captured.out

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
                    "left": ("AT", "GT", 2),
                    "right": ("GC", "GT", 4)
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
                    "left": ("AT", "GT", 2),
                    "right": ("GT", "GT", 0)
                }
            },
        }

        expected_result = {
            ('+', 'start', 1, 'left'): {
                ("AT", "GT"): {2: 2}
            },
            ('+', 'start', 1, 'right'): {
                ("GC", "GT"): {4: 1},
                ("GT", "GT"): {0: 1}
            }
        }
        self.assertEqual(
            log_manager.compute_closest_canonicals_dict(), expected_result)

    def test_errors_are_written_to_file(self):
        if os.path.exists(log_manager.alignment_error_log_filepath):
            os.remove(log_manager.alignment_error_log_filepath)
        log_manager.alignment_erros.append("test_error")
        log_manager.write_alignment_errors_to_file()
        self.assertTrue(os.path.exists(
            log_manager.alignment_error_log_filepath))

    def test_validation_for_indel_graphs_returns_true_if_enough_cases_and_closest_canonical_matches(self):
        key = ("+", "start", 1, "left")
        nucleotides = ("AG", "AC")
        cases = {2: 10}
        result = log_manager.validate_img_should_be_created_for_closest_canonical_dict_entry(
            key, nucleotides, cases)
        self.assertTrue(result)

    def test_validation_for_indel_graphs_returns_false_if_too_few_cases(self):
        key = ("+", "start", 1, "left")
        nucleotides = ("AT", "GT")
        cases = {2: 1}
        result = log_manager.validate_img_should_be_created_for_closest_canonical_dict_entry(
            key, nucleotides, cases)
        self.assertFalse(result)

    def test_validation_for_indel_graphs_returns_false_if_no_canonical_match(self):
        key = ("+", "start", 1, "left")
        nucleotides = ("AA", "GT")
        cases = {2: 10}
        result = log_manager.validate_img_should_be_created_for_closest_canonical_dict_entry(
            key, nucleotides, cases)
        self.assertFalse(result)

    def test_closest_canonical_generation_runs_correctly(self):
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
                    "left": ("AT", "GT", 2),
                    "right": ("GC", "GT", 4)
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
                    "left": ("AT", "GT", 2),
                    "right": ("GT", "GT", 0)
                }
            },
        }

        log_manager.generate_output_for_closest_canonicals(self.parser)

        captured = self.capsys.readouterr()
        assert "2" in captured.out

    def test_indel_graph_validator_returns_true_if_enough_cases(self):
        cases = {2: 11}
        self.parser.min_reads_for_graph = 10
        result = log_manager.validate_indel_grahp_should_be_created(
            self.parser, cases)
        self.assertTrue(result)

    def test_indel_graph_validator_returns_true_if_not_enough_cases(self):
        cases = {2: 9}
        self.parser.min_reads_for_graph = 10
        result = log_manager.validate_indel_grahp_should_be_created(
            self.parser, cases)
        self.assertFalse(result)

    def tearDown(self):
        if os.path.exists(log_manager.alignment_error_log_filepath):
            os.remove(log_manager.alignment_error_log_filepath)


class TestOffsetComputationFileManagement(TestCase):

    def setUp(self) -> None:
        if os.path.exists(OFFSET_LOG):
            os.remove(OFFSET_LOG)
        self.log_manager = LogManager()

    def test_a_test_file_is_initialized(self):
        self.log_manager.offset_results = {"test": {
            "reference_id": "test",
            "strand": "+",
            "class_code": "j",
            "offsets": [(1, 2), (3, 4)],
        }}
        self.log_manager.write_offset_results_to_file()
        self.assertTrue(os.path.exists(OFFSET_LOG))

    def test_initialized_log_file_only_has_header(self):
        self.log_manager.offset_results = {"test": {
            "reference_id": "test",
            "strand": "+",
            "class_code": "j",
            "offsets": [(1, 2), (3, 4)],
        }}
        self.log_manager.write_offset_results_to_file()
        with open(OFFSET_LOG, 'r') as f:
            lines = f.readlines()
            print(lines)
            self.assertEqual(len(lines), 2)

    def test_lines_are_written_to_log_file(self):
        self.log_manager.offset_results = {
            "transcript_1": {
                "reference_id": "ref_1",
                "strand": "+",
                "class_code": "j",
                "offsets": [(1, 2), (3, 4)],
            },
            "transcript_2": {
                "reference_id": "ref_2",
                "strand": "-",
                "class_code": "j",
                "offsets": [(1, 2), (3, 4)],
            },
            "transcript_3": {
                "reference_id": "ref_3",
                "strand": "+",
                "class_code": "x",
                "offsets": [(1, 2), (3, 4)],
            },
        }

        self.log_manager.write_offset_results_to_file()
        with open(OFFSET_LOG, 'r') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 4)

    def tearDown(self) -> None:
        if os.path.exists(OFFSET_LOG):
            os.remove(OFFSET_LOG)
