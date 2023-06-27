from services.bam_manager import BamManager

from services.offset_computation import execute_offset_computation
from services.argument_parser import init_argparser
from services.db_initializer import init_databases
from services.class_code_stats import ClassCodeStats
from services.extract_matching_cases import MatchingCasesExtractor
from services.fasta_extractor import FastaExtractor
from services.output_manager import default_output_manager as output_manager


def run_pipeline(parser_args):
    output_manager.output_heading()

    gffcompare_db, reference_db = init_databases(
        parser_args.gffcompare_gtf, parser_args.reference_gtf, parser_args.force)

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

    if matching_cases_dict and parser_args.reference_fasta:
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

    if matching_cases_dict and parser_args.reads_tsv and parser_args.reads_bam:
        bam_manager = BamManager(
            parser_args.reads_bam,
            parser_args.reads_tsv,
            matching_cases_dict,
            parser_args.extended_debug
        )
        bam_manager.execute(parser_args.window_size)
    output_manager.output_footer()
    output_manager.write_log_file()


def main():
    parser_args = init_argparser()
    run_pipeline(parser_args)


if __name__ == "__main__":
    main()
