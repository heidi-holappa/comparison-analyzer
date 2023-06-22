import sys
from unittest import TestCase

from compana import main


class TestCompana(TestCase):

    def setUp(self):
        self.execute = main
        pass

    def test_compana_runs_without_errors(self):
        sys.argv = [
            'compAna.py',
            '-j=src/tests/files/test-run.json'
        ]
        self.execute()
