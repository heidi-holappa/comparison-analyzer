import os
from unittest import TestCase

from services.bam_manager import BamManager
from config import TEMPORARY_DIR


class TestBamManager(TestCase):

    def setUp(self):
        self.bam_path = "bam_path"
        self.tsv_path = "tsv_path"
        self.matching_cases_dict = {}
        self.bam_manager = BamManager(
            self.bam_path, self.tsv_path, self.matching_cases_dict)

    def test_bam_manager_is_initialized_correctly(self):
        temporary_path = TEMPORARY_DIR
        self.assertEqual(self.bam_manager.bam_path, self.bam_path)
        self.assertEqual(self.bam_manager.tsv_path, self.tsv_path)
        self.assertEqual(self.bam_manager.matching_cases_dict,
                         self.matching_cases_dict)
        self.assertEqual(self.bam_manager.transcript_set, set())
        self.assertEqual(self.bam_manager.temporary_path, temporary_path)

    def test_temporary_directory_is_created(self):
        self.assertFalse(os.path.exists(self.bam_manager.temporary_path))
        self.bam_manager.create_temporary_path()
        self.assertTrue(os.path.exists(self.bam_manager.temporary_path))

    def test_temporary_directory_is_removed(self):
        self.bam_manager.create_temporary_path()
        self.assertTrue(os.path.exists(self.bam_manager.temporary_path))
        self.bam_manager.remove_temporary_path()
        self.assertFalse(os.path.exists(self.bam_manager.temporary_path))

    def tearDown(self) -> None:
        if os.path.exists(self.bam_manager.temporary_path):
            self.bam_manager.remove_temporary_path()
