from unittest import TestCase

from implementation_services.read_extractor import create_dict_of_transcripts_and_reads
from implementation_services.read_extractor import create_dict_of_reads_and_references
from tests.sample_file_management import default_test_file_manager


class TestImpReadExtractor(TestCase):

    def test_reads_are_correctly_extracted_from_file(self):
        transcripts = {'transcript1.chr1.nnic'}
        tsv_file = default_test_file_manager.tsv_file
        result = create_dict_of_transcripts_and_reads(tsv_file, transcripts)
        self.assertEqual(len(result), 1)

    def test_reads_and_references_are_correctly_created(self):
        intron_site_dict = {
            'intron_site1': {
                'transcript_id': 'transcript1.chr1.nnic'
            },
            'intron_site2': {
                'transcript_id': 'transcript3.chr1.nnic'
            }
        }
        dict_of_transcripts_and_reads = {
            'transcript1.chr1.nnic': {'read1'},
            'transcript2.chr1.nnic': {'read2'}
        }
        result = create_dict_of_reads_and_references(
            intron_site_dict, dict_of_transcripts_and_reads)
        self.assertEqual(len(result), 1)
