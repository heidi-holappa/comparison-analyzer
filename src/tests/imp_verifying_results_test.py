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
                'offset': 2
            },
            'key4': {
                'offset': 2
            }
        }

        intron_site_dict = {
            'key1': {
                'extracted_information': {
                    'right': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2)
                    },
                    'left': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2)
                    }
                },
            },
            'key2': {
                'extracted_information': {
                    'right': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2)
                    },
                    'left': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2)
                    }
                },
            },
            'key3': {
                'extracted_information': {
                    'right': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2)
                    },
                    'left': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2)
                    }
                },
            },
            'key4': {
                'extracted_information': {
                    'right': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2)
                    },
                    'left': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2)
                    }
                },
            },
            'key5': {
                'extracted_information': {
                    'right': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 1)
                    },
                    'left': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 1)
                    }
                },
            }
        }

        verify_results(intron_site_dict, matching_cases_dict)
        captured = self.capsys.readouterr()
        assert "True positives: {'left': {2: 1}, 'right': {2: 1}}" in captured.out
        assert "False positives: {'left': {}, 'right': {2: 1}}" in captured.out
