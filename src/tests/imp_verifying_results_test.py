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
                'offset': 2
            },
            'key3': {
                'offset': 4
            },
            'key4': {
                'offset': 4
            }
        }

        intron_site_dict = {
            'key1': {
                'extracted_information': {
                    'right': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    },
                    'left': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    }
                },
            },
            'key2': {
                'extracted_information': {
                    'right': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 4: 100},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3, 3, 3, 100, 100, 100, 100],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    },
                    'left': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    }
                },
            },
            'key3': {
                'extracted_information': {
                    'right': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    },
                    'left': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3, 4: 100},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 3),
                        'del_pos_distr': [1, 2, 3, 3, 3, 100, 100, 100, 100],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    }
                },
            },
            'key4': {
                'extracted_information': {
                    'right': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3, 4: 100},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 4),
                        'del_pos_distr': [1, 2, 3, 3, 3, 100, 100, 100, 100],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    },
                    'left': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    }
                },
            },
            'key5': {
                'extracted_information': {
                    'right': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3, 4: 100},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 1),
                        'del_pos_distr': [1, 2, 3, 3, 3, 100, 100, 100, 100],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    },
                    'left': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 1),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966
                    }
                },
            }
        }

        verify_results(intron_site_dict, matching_cases_dict)
        captured = self.capsys.readouterr()
        print(captured.out)
        assert "True positives: {'left': {4: 1}, 'right': {4: 1}, 'closest_canonical_matches': 1" in captured.out
        assert "False positives: {'left': {}, 'right': {4: 1}, 'closest_canonical_matches': 1" in captured.out
