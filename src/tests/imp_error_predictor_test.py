import sys

from unittest import TestCase
import pytest

from services.argument_parser import init_argparser

from implementation_services.error_predictor import execute_error_prediction
from implementation_services.error_predictor import make_prediction
from implementation_services.error_predictor import count_predicted_errors
from implementation_services.error_predictor import verify_sublist_largest_values_exists


class TestErrorPredictor(TestCase):

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

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_execute_prediction_gives_correct_output(self):
        intron_site_dict = {}
        execute_error_prediction(self.parser, intron_site_dict)
        captured = self.capsys.readouterr()
        assert "error prediction finished" in captured.out

    def test_prediction_maker_returns_if_findings_n_too_low(self):
        findings = {
            'insertions': {0: 0},
            'deletions': {0: 0},
            'closest_canonical': ('AC', 'AC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 0, 0, 0, 0],
        }

        expected_result = {
            'insertions': {0: 0},
            'deletions': {0: 0},
            'closest_canonical': ('AC', 'AC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 0, 0, 0, 0],
        }
        make_prediction(self.parser, findings, 'start', '+')
        self.assertEqual(findings, expected_result)

    def test_prediction_finds_error_for_exon_end(self):
        findings = {
            'insertions': {1: 100},
            'deletions': {4: 100},
            'closest_canonical': ('GC', 'GC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 80, 80, 80, 80],
            'most_common_del_pair': 'GT',
        }

        expected_result = {
            'error_detected': True,
        }
        make_prediction(self.parser, findings, 'end', '+')
        self.assertEqual(findings['error_detected'],
                         expected_result['error_detected'])

    def test_prediction_find_error_for_exon_start(self):
        findings = {
            'insertions': {1: 100},
            'deletions': {4: 100},
            'closest_canonical': ('AC', 'AC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 80, 80, 80, 80],
            'most_common_del_pair': 'AG',
        }

        expected_result = {
            'error_detected': True,
        }
        make_prediction(self.parser, findings, 'start', '+')
        self.assertEqual(findings['error_detected'],
                         expected_result['error_detected'])

    def test_avg_and_sd_are_calculated_correctly_for_insertions(self):
        findings = {
            'insertions': {1: 100, 2: 100, 3: 100, 4: 100, 5: 100},
            'deletions': {4: 100},
            'closest_canonical': ('AC', 'AC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 80, 80, 80, 80],
            'most_common_del_pair': 'AC',
        }
        expected_result = {
            'ins_avg': 3.0,
            'ins_sd': 1.4142135623730951,
            'del_avg': 4.0,
            'del_sd': 0.0,
        }
        make_prediction(self.parser, findings, 'start', '+')
        self.assertEqual(findings['ins_avg'], expected_result['ins_avg'])
        self.assertEqual(findings['ins_sd'], expected_result['ins_sd'])

    def test_only_one_error_is_counted_per_intron_site(self):
        intron_site_dict = {
            'key1': {
                'extracted_information': {
                    'left': {
                        'error_detected': True
                    },
                    'right': {
                        'error_detected': True
                    }
                }
            },
        }
        count_predicted_errors(intron_site_dict)
        captured = self.capsys.readouterr()
        assert "Predicted errors: 1" in captured.out

    def test_if_no_sublist_of_largest_values_exists_return_false(self):
        l = [10, 10, 0, 0, 10, 10, 0, 0, 0]
        n = 3
        result = verify_sublist_largest_values_exists(l, n)
        self.assertEqual(result, False)

    def test_if_canonicals_strategy_and_canonicals_do_not_match_no_error_is_detected(self):
        self.parser.no_canonicals = False
        findings = {
            'insertions': {1: 100, 2: 100, 3: 100, 4: 100, 5: 100},
            'deletions': {4: 100},
            'closest_canonical': ('AC', 'AC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 80, 80, 80, 80],
            'most_common_del_pair': 'GG',
        }

        make_prediction(self.parser, findings,
                        location_type='start', strand='+')
        self.assertEqual(False, findings['error_detected'])

    def test_if_not_accepted_offset_case_do_nothing(self):
        findings = {
            'insertions': {1: 100, 2: 100, 3: 100, 4: 100, 5: 100},
            'deletions': {1: 100},
            'closest_canonical': ('AC', 'AC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 80, 80, 80, 80],
            'most_common_del_pair': 'GG',
        }

        make_prediction(self.parser, findings,
                        location_type='start', strand='+')
        self.assertEqual(False, findings['error_detected'])

    def test_very_conservative_strategy_returns_true_if_conditions_are_met(self):
        # For conditions see documentation
        findings = {
            'insertions': {0: 100},
            'deletions': {4: 100},
            'closest_canonical': ('AC', 'AC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 80, 80, 80, 80],
            'most_common_del_pair': 'AG',
        }
        self.parser.no_canonicals = False
        make_prediction(self.parser, findings,
                        location_type='start', strand='+')
        self.assertEqual(True, findings['error_detected'])
