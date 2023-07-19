import sys

from services.bam_manager import BamManager

from services.offset_computation import execute_offset_computation
from services.argument_parser import init_argparser
from services.db_initializer import init_databases
from services.class_code_stats import ClassCodeStats
from services.extract_matching_cases import MatchingCasesExtractor
from services.fasta_extractor import FastaExtractor
from services.output_manager import default_output_manager as output_manager
from services.log_manager import default_log_manager as log_manager

from implementation_services.isoquant_db_init import init_isoquant_db
from implementation_services.case_extractor import CaseExtractor
from implementation_services.read_extractor import create_dict_of_transcripts_and_reads, create_dict_of_reads_and_references
from implementation_services.indel_computer import execute_indel_computation


def run_pipeline(parser_args):
    output_manager.output_heading()

    gffcompare_db, reference_db = init_databases(
        parser_args.gffcompare_gtf,
        parser_args.reference_gtf,
        parser_args.force)

    if parser_args.stats:
        class_code_stats = ClassCodeStats(gffcompare_db)
        class_code_stats.compute_class_code_stats()

    offset_results = {}
    if parser_args.class_code:
        offset_results = execute_offset_computation(
            parser_args.class_code, gffcompare_db, reference_db, parser_args.extended_debug)

    extractor = MatchingCasesExtractor(
        offset_results,
        parser_args.offset,
        reference_db)
    matching_cases_dict = extractor.extract_candidates_matching_selected_offset()

    if parser_args.reference_fasta:
        fasta_config = {
            "fasta_path": parser_args.reference_fasta,
            "offset": parser_args.offset,
            "gffcompare_gtf": parser_args.gffcompare_gtf,
            "reference_gtf": parser_args.reference_gtf,
            "class_codes": parser_args.class_code,
            "matching_cases_dict": matching_cases_dict,
            "window_size": parser_args.window_size
        }
        reference_fasta_extractor = FastaExtractor(fasta_config)
        reference_fasta_extractor.execute_fasta_extraction()

    if parser_args.reads_tsv and parser_args.reads_bam:
        bam_manager = BamManager(
            parser_args.reads_bam,
            parser_args.reads_tsv,
            matching_cases_dict)
        bam_manager.execute(int(parser_args.window_size))

    log_manager.execute_log_file_creation(matching_cases_dict, parser_args)

    if not parser_args.isoquant_gtf:
        sys.exit(1)

    isoquant_db = init_isoquant_db(parser_args.isoquant_gtf)

    # Extract all intron site locations from isoquant-db and create a set of transcripts
    # input: isoquant-db
    # output: intron site dictionary

    case_extractor = CaseExtractor()
    intron_site_dict = case_extractor.extract_intron_site_locations(
        isoquant_db)

    transcript_set = case_extractor.extract_transcripts(isoquant_db)

    # extract reads and references
    # input model_reads.tsv, intron site dictionary
    # output: reads and references dictionary

    transcripts_and_reads = create_dict_of_transcripts_and_reads(
        parser_args.reads_tsv, transcript_set)

    reads_and_references = create_dict_of_reads_and_references(
        intron_site_dict, transcripts_and_reads)

    # compute indels
    # input: intorn site dictionary, reads and references dictionary, bam file
    # output: updated intron site dictionary

    execute_indel_computation(
        parser_args.reads_bam,
        intron_site_dict,
        reads_and_references,
        parser_args.window_size)

    # compute closest canonicals for interesting cases
    # input: intron site dictionary, reference fasta file
    # output: updated intron site dictionary

    # Predict possible mistakes based on indels and closest canonicals
    # Input: intron site dictionary
    # Output: updated transcript model

    # verify results
    # input: intron site dictionary, matching cases dict
    # output: verification results: misses, hits, errors


def main():
    parser_args = init_argparser()
    run_pipeline(parser_args)


if __name__ == "__main__":
    main()
