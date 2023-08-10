import sys

from unittest import TestCase
import pytest

from services.argument_parser import init_argparser

from implementation_services.verify_results import verify_results


class TestVerifyResults(TestCase):

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
        ]
        self.parser = init_argparser()

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
                'strand': '+',
                'extracted_information': {
                    'right': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    },
                    'left': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    }
                },
            },
            'key2': {
                'strand': '+',
                'extracted_information': {
                    'right': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 4: 100},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3, 3, 3, 100, 100, 100, 100],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    },
                    'left': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    }
                },
            },
            'key3': {
                'strand': '+',
                'extracted_information': {
                    'right': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    },
                    'left': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3, 4: 100},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 3),
                        'del_pos_distr': [1, 2, 3, 3, 3, 100, 100, 100, 100],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    }
                },
            },
            'key4': {
                'strand': '+',
                'extracted_information': {
                    'right': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3, 4: 100},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 4),
                        'del_pos_distr': [1, 2, 3, 3, 3, 100, 100, 100, 100],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    },
                    'left': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 2),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    }
                },
            },
            'key5': {
                'strand': '+',
                'extracted_information': {
                    'right': {
                        'error_detected': True,
                        'deletions': {0: 1, 1: 2, 2: 3, 4: 100},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 1),
                        'del_pos_distr': [1, 2, 3, 3, 3, 100, 100, 100, 100],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    },
                    'left': {
                        'error_detected': False,
                        'deletions': {0: 1, 1: 2, 2: 3},
                        'insertions': {0: 1, 1: 2, 2: 3},
                        'closest_canonical': ('CG', 'GT', 1),
                        'del_pos_distr': [1, 2, 3],
                        'del_avg': 2.0033333,
                        'del_sd': 0.8164966,
                        'most_common_del_pair': 'GT',
                    }
                },
            }
        }

        verify_results(self.parser, intron_site_dict, matching_cases_dict)
        captured = self.capsys.readouterr()
        print(captured.out)
        assert "True positives: {'left': {4: 1}, 'right': {4: 1}, 'closest_canonical_matches': 1" in captured.out
        assert "False positives: {'left': {}, 'right': {4: 1}, 'closest_canonical_matches': 1" in captured.out
