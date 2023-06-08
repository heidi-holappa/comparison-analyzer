import os
from unittest import TestCase
from unittest.mock import MagicMock, patch
import pytest

from services.db_initializer import init_databases
from services import db_initializer
from config import TEST_FILE_DIR


class TestCompana(TestCase):

    def setUp(self):
        self.gff_gtf = os.path.join(TEST_FILE_DIR, "gffcompare.annotated.gtf")
        self.ref_gtf = os.path.join(
            TEST_FILE_DIR, "gencode.vM26.basic.annotation.extract.gtf")

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    @patch('services.db_initializer.init_databases')
    def test_init_databases_is_called(self, mock_init_databases):
        mock_init_databases.return_value = MagicMock(), MagicMock()
        db_initializer.init_databases(MagicMock(), MagicMock(), False)
        mock_init_databases.assert_called_once()

    def test_databases_are_initialized(self):

        init_databases(self.gff_gtf, self.ref_gtf, True)
        self.assertTrue(os.path.exists(
            os.path.join(TEST_FILE_DIR, "gffcompare.annotated-ca.db")
        )
        )
        self.assertTrue(os.path.exists(
            os.path.join(
                TEST_FILE_DIR, "gencode.vM26.basic.annotation.extract-ca.db")
        )
        )

    def test_if_databases_exists_and_force_not_used_do_not_initialize(self):
        init_databases(self.gff_gtf, self.ref_gtf, True)
        init_databases(self.gff_gtf, self.ref_gtf, False)
        captured = self.capsys.readouterr()
        assert "using existing db" in captured.out

    def tearDown(self):
        if os.path.exists(os.path.join(TEST_FILE_DIR, "gffcompare.annotated-ca.db")):
            os.remove(os.path.join(TEST_FILE_DIR,
                      "gffcompare.annotated-ca.db"))
        if os.path.exists(os.path.join(TEST_FILE_DIR, "gencode.vM26.basic.annotation.extract-ca.db")):
            os.remove(os.path.join(
                TEST_FILE_DIR, "gencode.vM26.basic.annotation.extract-ca.db"))
