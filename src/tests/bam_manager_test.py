import os
from unittest import TestCase

from services.bam_manager import BamManager
from tests.sample_file_management import default_test_file_manager as file_manager


class TestBamManagerInit(TestCase):

    def setUp(self):
        self.bam_path = file_manager.bam_file
        self.tsv_path = file_manager.tsv_file
        self.matching_cases_dict = {
            'transcript1.chr1.nnic.exon_4.start': {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '+',
                'location_type': 'start',
                'exon_number': 3,
                'location': 210,
            },
            'transcript1.chr1.nnic.exon_4.end': {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '+',
                'location_type': 'end',
                'exon_number': 4,
                'location': 45,
            },
            'transcript3.chr1.nnic.exon_2.start': {
                'transcript_id': 'transcript3.chr1.nnic',
                'strand': '-',
                'location_type': 'start',
                'exon_number': 1,
                'location': 20,
            }
        }
        self.bam_manager = BamManager(
            self.bam_path, self.tsv_path, self.matching_cases_dict)

    def test_bam_manager_is_initialized_correctly(self):
        self.assertEqual(self.bam_manager.bam_path, self.bam_path)
        self.assertEqual(self.bam_manager.tsv_path, self.tsv_path)
        self.assertEqual(self.bam_manager.matching_cases_dict,
                         self.matching_cases_dict)
        self.assertEqual(self.bam_manager.transcript_set, set(
            ["transcript1.chr1.nnic", "transcript3.chr1.nnic"]))


class TestBamManagerExecution(TestCase):

    def setUp(self):
        self.bam_path = file_manager.bam_file
        self.tsv_path = file_manager.tsv_file
        self.extended_debugging = True
        self.matching_cases_dict = {
            'transcript1.chr1.nnic.exon_4.start': {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '+',
                'location_type': 'start',
                'exon_number': 3,
                'location': 210,
            },
        }
        self.bam_manager = BamManager(
            self.bam_path, self.tsv_path, self.matching_cases_dict)

    # TODO: activate test

    def test_bam_manager_execute_runs_without_errors(self):
        window_size = 10
        self.bam_manager.execute(window_size)


class TestDebugLogWriting(TestCase):

    def setUp(self):
        self.bam_manager = BamManager("", "", {})
        if os.path.exists(self.bam_manager.debug_log_path_reads_and_locations):
            os.remove(self.bam_manager.debug_log_path_reads_and_locations)
        if os.path.exists(self.bam_manager.debug_log_path_transcripts_and_reads):
            os.remove(self.bam_manager.debug_log_path_transcripts_and_reads)

    def test_debug_log_is_written(self):
        self.bam_manager.write_debug_logs({"test": "test"}, {"test": "test"})
        self.assertTrue(os.path.exists(
            self.bam_manager.debug_log_path_transcripts_and_reads))
        self.assertTrue(os.path.exists(
            self.bam_manager.debug_log_path_reads_and_locations))

    def tearDown(self):
        if os.path.exists(self.bam_manager.debug_log_path_reads_and_locations):
            os.remove(self.bam_manager.debug_log_path_reads_and_locations)
        if os.path.exists(self.bam_manager.debug_log_path_transcripts_and_reads):
            os.remove(self.bam_manager.debug_log_path_transcripts_and_reads)
