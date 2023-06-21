import os
from unittest import TestCase
from unittest.mock import MagicMock, patch
import pytest

from services.db_initializer import init_databases
from tests.sample_file_management import default_test_file_manager
from services import db_initializer
from config import TEST_FILE_DIR


class TestCompana(TestCase):

    def setUp(self):
        self.test_file_manager = default_test_file_manager
        self.test_file_manager.initialize_test_files()

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    @patch('services.db_initializer.init_databases')
    def test_init_databases_is_called(self, mock_init_databases):
        mock_init_databases.return_value = MagicMock(), MagicMock()
        db_initializer.init_databases(MagicMock(), MagicMock(), False)
        mock_init_databases.assert_called_once()

    def test_databases_are_initialized(self):
        init_databases(self.test_file_manager.gffcompare_gtf,
                       self.test_file_manager.reference_gtf,
                       True)
        self.assertTrue(self.test_file_manager.databases_exist())

    def test_if_databases_exists_and_force_not_used_do_not_initialize(self):
        init_databases(self.test_file_manager.gffcompare_gtf,
                       self.test_file_manager.reference_gtf,
                       True)
        init_databases(self.test_file_manager.gffcompare_gtf,
                       self.test_file_manager.reference_gtf,
                       False)
        captured = self.capsys.readouterr()
        assert "exists" in captured.out

    def tearDown(self):
        self.test_file_manager.remove_test_files()
