# Consider how to best test the main function
# The issue is that now this test runs the whole pipeline
# and 'covers' most of the code. However, it is not a unit test.
# It is more of an integration test.
# Consider mocking the functions that are called in main

import sys
from unittest import TestCase

from compana import main


class TestCompana(TestCase):

    def setUp(self):
        self.execute = main
        pass

    def test_compana_runs_without_errors_with_all_args(self):
        sys.argv = [
            'compAna.py',
            '-j=src/tests/files/test-run.json'
        ]
        self.execute()

    def test_compana_runs_without_errors_with_limited_args(self):
        sys.argv = [
            'compAna.py',
            '-j=src/tests/files/test-run-fewer-args.json'
        ]
        self.execute()
