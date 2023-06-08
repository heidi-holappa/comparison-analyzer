import os
from unittest import TestCase
import pytest

from services.class_code_stats import compute_class_code_stats
from services.db_initializer import init_databases
from services import db_initializer
from config import TEST_FILE_DIR


class TestClassCodeStats(TestCase):

    def setUp(self):
        self.gff_gtf = os.path.join(TEST_FILE_DIR, "gffcompare.annotated.gtf")
        self.ref_gtf = os.path.join(
            TEST_FILE_DIR, "gencode.vM26.basic.annotation.extract.gtf")
        self.gff_db, self.ref_db = init_databases(
            self.gff_gtf, self.ref_gtf, True)

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_class_code_stats_outputs_class_code_statistics_of_results(self):
        results = compute_class_code_stats(self.gff_db)
        captured = self.capsys.readouterr()
        assert "class code" in captured.out

    def tearDown(self):
        if os.path.exists(os.path.join(TEST_FILE_DIR, "gffcompare.annotated-ca.db")):
            os.remove(os.path.join(TEST_FILE_DIR,
                      "gffcompare.annotated-ca.db"))
        if os.path.exists(os.path.join(TEST_FILE_DIR, "gencode.vM26.basic.annotation.extract-ca.db")):
            os.remove(os.path.join(
                TEST_FILE_DIR, "gencode.vM26.basic.annotation.extract-ca.db"))
