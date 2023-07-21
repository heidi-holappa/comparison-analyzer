from unittest import TestCase
import pytest

from implementation_services.verify_results import verify_results


class TestVerifyResults(TestCase):

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_results_are_correctly_output_to_stdout(self):
        matching_cases_dict = {
            'key1': {
                'offset': 0
            },
            'key2': {
                'offset': 1
            },
            'key3': {
                'offset': 0
            },
            'key4': {
                'offset': 2
            }
        }

        intron_site_dict = {
            'key1': {
                'extracted_information': {
                    'right': {
                        'error_detected': False
                    },
                    'left': {
                        'error_detected': False
                    }
                },
            },
            'key2': {
                'extracted_information': {
                    'right': {
                        'error_detected': True
                    },
                    'left': {
                        'error_detected': False
                    }
                },
            },
            'key3': {
                'extracted_information': {
                    'right': {
                        'error_detected': False
                    },
                    'left': {
                        'error_detected': True
                    }
                },
            },
            'key4': {
                'extracted_information': {
                    'right': {
                        'error_detected': True
                    },
                    'left': {
                        'error_detected': True
                    }
                },
            },
            'key5': {
                'extracted_information': {
                    'right': {
                        'error_detected': True
                    },
                    'left': {
                        'error_detected': True
                    }
                },
            }
        }

        verify_results(intron_site_dict, matching_cases_dict)
        captured = self.capsys.readouterr()
        assert "True positives: 2" in captured.out
        assert "False positives: 1" in captured.out
