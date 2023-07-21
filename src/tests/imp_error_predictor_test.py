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
            'closest_canonical': ('', '', 0),
            'error_detected': False
        }

        expected_result = {
            'insertions': {0: 0},
            'deletions': {0: 0},
            'closest_canonical': ('', '', 0),
            'error_detected': False
        }
        make_prediction(findings)
        self.assertEqual(findings, expected_result)

    def test_prediction_maker_does_not_give_key_error_if_no_zeroes(self):
        findings = {
            'insertions': {1: 100},
            'deletions': {1: 100},
            'closest_canonical': ('', '', 4),
            'error_detected': False
        }

        expected_result = {
            'insertions': {1: 100},
            'deletions': {1: 100},
            'closest_canonical': ('', '', 4),
            'error_detected': True
        }
        make_prediction(findings)
        self.assertEqual(findings, expected_result)

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
