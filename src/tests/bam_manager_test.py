import os
from unittest import TestCase

from services.bam_manager import BamManager
from tests.sample_file_management import default_test_file_manager as file_manager
from config import TEMPORARY_DIR


class TestBamManagerInit(TestCase):

    def setUp(self):
        self.bam_path = file_manager.bam_file
        self.tsv_path = file_manager.tsv_file
        self.matching_cases_dict = {
            ("transcript_1", "case_1"): "case_1",
            ("transcript_2", "case_2"): "case_2",
            ("transcript_3", "case_3"): "case_3"
        }
        self.bam_manager = BamManager(
            self.bam_path, self.tsv_path, self.matching_cases_dict)

    def test_bam_manager_is_initialized_correctly(self):
        self.assertEqual(self.bam_manager.bam_path, self.bam_path)
        self.assertEqual(self.bam_manager.tsv_path, self.tsv_path)
        self.assertEqual(self.bam_manager.matching_cases_dict,
                         self.matching_cases_dict)
        self.assertEqual(self.bam_manager.transcript_set, set(
            ["transcript_1", "transcript_2", "transcript_3"]))

    def test_temporary_directory_is_created(self):
        self.assertFalse(os.path.exists(TEMPORARY_DIR))
        self.bam_manager.create_temporary_path()
        self.assertTrue(os.path.exists(TEMPORARY_DIR))

    def test_temporary_directory_is_removed(self):
        self.bam_manager.create_temporary_path()
        self.assertTrue(os.path.exists(TEMPORARY_DIR))
        self.bam_manager.remove_temporary_path()
        self.assertFalse(os.path.exists(TEMPORARY_DIR))

    def tearDown(self) -> None:
        if os.path.exists(TEMPORARY_DIR):
            self.bam_manager.remove_temporary_path()


class TestBamManagerExecution(TestCase):

    def setUp(self):
        self.bam_path = file_manager.bam_file
        self.tsv_path = file_manager.tsv_file
        self.matching_cases_dict = {
            ("transcript1.chr1.nnic", "case_1"): "case_1"
        }
        self.bam_manager = BamManager(
            self.bam_path, self.tsv_path, self.matching_cases_dict)
        self.bam_manager.create_temporary_path()

    def test_bam_manager_execute_runs_without_errors(self):
        self.bam_manager.execute()

    def tearDown(self) -> None:
        if os.path.exists(TEMPORARY_DIR):
            self.bam_manager.remove_temporary_path()
