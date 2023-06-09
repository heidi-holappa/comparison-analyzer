from services.bam_manager import BamManager

from services.offset_computation import execute_offset_computation
from services.argument_parser import init_argparser
from services.db_initializer import init_databases
from services.class_code_stats import ClassCodeStats
from services.fasta_extractor import execute_fasta_extraction
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
            parser_args.class_code, gffcompare_db, reference_db)

    matching_cases_dict = {}

    if parser_args.reference_fasta:
        matching_cases_dict = execute_fasta_extraction(
            parser_args, offset_results, reference_db)

    if parser_args.reads_tsv and parser_args.reads_bam:
        bam_manager = BamManager(
            parser_args.reads_bam,
            parser_args.reads_tsv,
            matching_cases_dict
        )
        bam_manager.execute()
    output_manager.output_footer()


def main():
    parser_args = init_argparser()
    run_pipeline(parser_args)


if __name__ == "__main__":
    main()
