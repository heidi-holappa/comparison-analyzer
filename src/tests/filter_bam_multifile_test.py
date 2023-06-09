import os
import sys
import shutil
from unittest import TestCase
import pytest

from services.filter_bam_multifile import create_read_dict
from services.filter_bam_multifile import create_output_filename_dict_cli
from services.filter_bam_multifile import create_output_filename_dict
from services.filter_bam_multifile import filter_reads
from services.filter_bam_multifile import init_parser
from tests.sample_file_management import default_test_file_manager as file_manager

from config import TEMPORARY_DIR


class TestReadDict(TestCase):
    def setUp(self):
        if not os.path.exists(TEMPORARY_DIR):
            os.makedirs(TEMPORARY_DIR)
        self.read_file = os.path.join(TEMPORARY_DIR, "read_list.tsv")
        with open(self.read_file, "w", encoding="UTF-8") as file:
            file.write("read1\ttranscript1\n")
            file.write("read1\ttranscript2\n")
            file.write("read2\ttranscript2\n")
            file.write("read3\ttranscript1\n")
            file.write("read4\ttranscript3\n")
        self.output_filename_dict = {
            "transcript1": "transcript1.bam",
            "transcript2": "transcript2.bam"
        }

    def test_read_dict_has_the_correct_number_of_keys(self):
        read_dict = create_read_dict(self.output_filename_dict, self.read_file)
        self.assertEqual(len(read_dict), 3)

    def test_read_dict_read1_has_two_transcripts(self):
        read_dict = create_read_dict(self.output_filename_dict, self.read_file)
        self.assertEqual(len(read_dict["read1"]), 2)

    def tearDown(self) -> None:
        if os.path.exists(TEMPORARY_DIR):
            shutil.rmtree(TEMPORARY_DIR)


class TestCreateFilenameDictionary(TestCase):

    def setUp(self):
        if not os.path.exists(TEMPORARY_DIR):
            os.makedirs(TEMPORARY_DIR)
        self.transcript_list_file = os.path.join(
            TEMPORARY_DIR, "transcript_list.txt")
        with open(self.transcript_list_file, "w", encoding="UTF-8") as file:
            file.write("transcript1\n")
            file.write("transcript2\n")
            file.write("transcript3\n")

    def test_output_filename_dict_has_correct_number_of_keys(self):
        bam_file = os.path.join(TEMPORARY_DIR, "test.bam")
        transcript_set = {"transcript1", "transcript2"}
        output_filename_dict = create_output_filename_dict(
            bam_file, transcript_set, TEMPORARY_DIR)
        self.assertEqual(len(output_filename_dict), 2)

    def test_cli_version_output_filename_dict_has_correct_number_of_keys(self):
        bam_file = os.path.join(TEMPORARY_DIR, "test.bam")
        output_filename_dict = create_output_filename_dict_cli(
            bam_file,
            self.transcript_list_file,
            suffix="test"
        )
        self.assertEqual(len(output_filename_dict), 3)

    def tearDown(self) -> None:
        if os.path.exists(TEMPORARY_DIR):
            shutil.rmtree(TEMPORARY_DIR)


class TestFilterBamFile(TestCase):

    def setUp(self):
        if not os.path.exists(TEMPORARY_DIR):
            os.makedirs(TEMPORARY_DIR)
        self.file_manager = file_manager
        self.bam_file = self.file_manager.bam_file

    def test_filter_bam_file_produces_output_bam_file(self):
        output_filename = os.path.join(TEMPORARY_DIR, "transcript1.bam")
        output_filename_dict = {
            "transcript925.ch1.nnic": output_filename
        }
        dict_of_reads = {
            "ENSMUST00000208994_962_aligned_8455583_F_112_261_0": ["transcript925.ch1.nnic"],
        }
        filter_reads(self.bam_file, output_filename_dict, dict_of_reads)
        self.assertTrue(os.path.exists(output_filename))

    def tearDown(self) -> None:
        if os.path.exists(TEMPORARY_DIR):
            shutil.rmtree(TEMPORARY_DIR)


class TestFilterBamParser(TestCase):

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_parser_reads_short_arguments_correctly(self):
        sys.argv = [
            'filter_bam_multifile.py',
            '-i', 'input.bam',
            '-m', 'model.tsv',
            '-t', 'transcript_list.txt',
            '-s', 'suffix',

        ]
        parser = init_parser()
        parser = parser.parse_args()
        assert parser.input == 'input.bam'
        assert parser.model_reads_tsv == 'model.tsv'
        assert parser.transcript_list == 'transcript_list.txt'
        assert parser.suffix == 'suffix'
