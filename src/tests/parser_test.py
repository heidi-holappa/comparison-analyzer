import unittest
import pytest
import sys
import json

from services.argument_parser import init_argparser


class TestParser(unittest.TestCase):

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_parser_reads_short_arguments_correctly(self):
        sys.argv = [
            'compAna.py',
            '-g=gffcompare.gtf',
            '-r=reference.gtf',
            '-b=file.bam',
            '-a=file.fasta',
            '-t=file.tsv',
            '-o=10',
            '-f',
            '-s',
            '-c=j k',

        ]
        parser = init_argparser()
        assert parser.gffcompare_gtf == 'gffcompare.gtf'
        assert parser.reference_gtf == 'reference.gtf'
        assert parser.reads_bam == 'file.bam'
        assert parser.reads_tsv == 'file.tsv'
        assert parser.offset == '10'
        assert parser.force
        assert parser.stats
        assert parser.class_code == ['j k']

    def test_parser_reads_long_arguments_correctly(self):
        sys.argv = [
            'compAna.py',
            '--gffcompare_gtf=gffcompare.gtf',
            '--reference_gtf=reference.gtf',
            '--reads_bam=file.bam',
            '--reference_fasta=file.fasta',
            '--reads_tsv=file.tsv',
            '--offset=10',
            '--force',
            '--stats',
            '--class-code=j k',

        ]
        parser = init_argparser()
        assert parser.gffcompare_gtf == 'gffcompare.gtf'
        assert parser.reference_gtf == 'reference.gtf'
        assert parser.reads_bam == 'file.bam'
        assert parser.reads_tsv == 'file.tsv'
        assert parser.offset == '10'
        assert parser.force
        assert parser.stats
        assert parser.class_code == ['j k']

    def test_parser_reads_json_arguments_correctly(self):
        sys.argv = [
            'compAna.py',
            '--json=src/tests/test.json'
        ]
        parser = init_argparser()
        with open('src/tests/test.json', encoding="UTF-8") as json_file:
            json_file = json.load(json_file)
            print(json_file)
        assert parser.gffcompare_gtf == 'gffcompare.gtf'
        assert parser.reference_gtf == 'reference.gtf'
        assert parser.reads_bam == 'file.bam'
        assert parser.reads_tsv == 'file.tsv'
        assert parser.offset == 10
        assert not parser.force
        assert parser.stats
        assert parser.class_code == 'j c k'

    def test_missing_args_returns_help(self):
        sys.argv = [
            'compAna.py',
        ]
        with self.assertRaises(SystemExit):
            init_argparser()
