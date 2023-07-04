import os
from unittest import TestCase

from services.bam_manager import BamManager
from tests.sample_file_management import default_test_file_manager as file_manager
from config import TEST_FILE_DIR


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


class TestBamManagerExecution(TestCase):

    def setUp(self):
        self.bam_path = file_manager.bam_file
        self.tsv_path = file_manager.tsv_file
        self.extended_debugging = True
        self.matching_cases_dict = {
            'transcript1.chr1.nnic.exon_4.start': {
                'transcript_id': 'transcript1.chr1.nnic',
                'strand': '+',
                'offset': 4,
                'location_type': 'start',
                'exon_number': 3,
                'location': 210,
            },
        }
        self.bam_manager = BamManager(
            self.bam_path, self.tsv_path, self.matching_cases_dict)

    def test_bam_manager_execute_runs_without_errors(self):
        window_size = 10
        self.bam_manager.execute(window_size)


class TestBamReader(TestCase):

    def setUp(self):
        self.bam_manager = BamManager("", "", {})
        test_bam_file = os.path.join(
            TEST_FILE_DIR, "Mouse.ONT.R9.4.sim.RE.no_gtf.transcript925.ch1.nnic.bam")
        self.bam_manager.initialize_file(test_bam_file)

    def test_bam_reader_runs_without_errors(self):
        reads_and_references = {
            "ENSMUST00000208994_1011_aligned_5112815_F_38_212_77": {'transcript1'}
        }
        self.bam_manager.matching_cases_dict = {
            'transcript1': {
                'location': 5112815,
                'location_type': 'start',
                'strand': '+',
                'offset': 4
            }

        }
        self.bam_manager.process_bam_file(reads_and_references)
