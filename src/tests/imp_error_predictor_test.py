from unittest import TestCase
import pytest

from implementation_services.error_predictor import execute_error_prediction
from implementation_services.error_predictor import make_prediction
from implementation_services.error_predictor import count_predicted_errors


class TestErrorPredictor(TestCase):

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_execute_prediction_gives_correct_output(self):
        intron_site_dict = {}
        execute_error_prediction(intron_site_dict)
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
        make_prediction(findings, 'start', '+')
        self.assertEqual(findings, expected_result)

    def test_prediction_finds_error_for_exon_end(self):
        findings = {
            'insertions': {1: 100},
            'deletions': {4: 100},
            'closest_canonical': ('GC', 'GC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 80, 80, 80, 80],
            'most_common_del_pair': (4, 'GT'),
        }

        expected_result = {
            'error_detected': True,
        }
        make_prediction(findings, 'end', '+')
        self.assertEqual(findings['error_detected'],
                         expected_result['error_detected'])

    def test_prediction_find_error_for_exon_start(self):
        findings = {
            'insertions': {1: 100},
            'deletions': {4: 100},
            'closest_canonical': ('AC', 'AC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 80, 80, 80, 80],
            'most_common_del_pair': (4, 'AG'),
        }

        expected_result = {
            'error_detected': True,
        }
        make_prediction(findings, 'start', '+')
        self.assertEqual(findings['error_detected'],
                         expected_result['error_detected'])

    def test_avg_and_sd_are_calculated_correctly_for_insertions(self):
        findings = {
            'insertions': {1: 100, 2: 100, 3: 100, 4: 100, 5: 100},
            'deletions': {4: 100},
            'closest_canonical': ('AC', 'AC', 4),
            'error_detected': False,
            'del_pos_distr': [0, 0, 0, 0, 80, 80, 80, 80],
            'most_common_del_pair': (4, 'AC'),
        }
        expected_result = {
            'ins_avg': 3.0,
            'ins_sd': 1.4142135623730951,
            'del_avg': 4.0,
            'del_sd': 0.0,
        }
        make_prediction(findings, 'start', '+')
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
