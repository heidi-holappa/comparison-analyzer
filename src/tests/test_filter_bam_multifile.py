import os
import shutil
from unittest import TestCase
from services.filter_bam_multifile import create_read_dict
from services.filter_bam_multifile import create_output_filename_dict_cli
from services.filter_bam_multifile import create_output_filename_dict

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
        print(output_filename_dict)
        self.assertEqual(len(output_filename_dict), 2)

    def test_cli_version_output_filename_dict_has_correct_number_of_keys(self):
        bam_file = os.path.join(TEMPORARY_DIR, "test.bam")
        output_filename_dict = create_output_filename_dict_cli(
            bam_file,
            self.transcript_list_file,
            suffix="test"
        )
        print(output_filename_dict)
        self.assertEqual(len(output_filename_dict), 3)

    def tearDown(self) -> None:
        if os.path.exists(TEMPORARY_DIR):
            shutil.rmtree(TEMPORARY_DIR)
