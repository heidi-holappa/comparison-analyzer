from services.bam_manager import BamManager

from services.offset_computation import execute_offset_computation
from services.argument_parser import init_argparser
from services.db_initializer import init_databases
from services.class_code_stats import compute_class_code_stats
from services.fasta_extractor import execute_fasta_extraction

parser = init_argparser()

print("=========================================")
print("compAna: a tool for comparing annotations")
print("=========================================")
print("===========DATABASE MANAGEMENT===========\n")

gffcompare_db, reference_db = init_databases(parser)


if parser.stats:
    compute_class_code_stats(gffcompare_db)

offset_results = {}
if parser.class_code:
    offset_results = execute_offset_computation(
        parser, gffcompare_db, reference_db)

matching_cases_dict = {}

if parser.reference_fasta:
    matching_cases_dict = execute_fasta_extraction(
        parser, offset_results, reference_db)

if parser.reads_tsv and parser.reads_bam:
    print("Fetching reads from BAM-file. This might take some time.")
    transcripts = set()
    for row in matching_cases_dict:
        transcripts.add(row[0])
    bam_manager = BamManager(
        parser.reads_bam, parser.reads_tsv, transcripts)
    bam_manager.execute()

print("================= END ===================")
