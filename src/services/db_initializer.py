import os
import gffutils
from services.output_manager import default_output_manager as output_manager


def init_databases(gffcompare_gtf: str, reference_gtf: str, force=False) -> tuple:
    gtf_paths = {
        'gffcompare': gffcompare_gtf,
        'reference': reference_gtf
    }

    db_paths = {
        'gffcompare': gffcompare_gtf[:-4] + '-ca.db',
        'reference': reference_gtf[:-4] + '-ca.db'
    }

    output_manager.output_line(
        "FILE INFORMATION", is_title=True)
    output_manager.output_line(
        f"Gffcompare GTF-file: {os.path.basename(gffcompare_gtf)}", is_info=True)
    output_manager.output_line(
        f"Reference GTF-file: {os.path.basename(reference_gtf)}", is_info=True)

    for key, value in db_paths.items():
        db_exists = os.path.exists(f'{value}')

        if not force and db_exists:
            output_manager.output_line(
                f"{key}: using existing db file. Use -f to force overwrite existing db-files.",
                is_info=True)
        else:
            output_manager.output_line(
                f'{key}: creating database... this might take a while.', is_info=True)
            gffutils.create_db(
                gtf_paths[key],
                dbfn=f'{value}',
                force=True,
                keep_order=True,
                merge_strategy='merge',
                sort_attribute_values=True,
                disable_infer_genes=True,
                disable_infer_transcripts=True
            )
            output_manager.output_line(
                f"{key}: database created successfully!", is_info=True)

    gffcompare_db = gffutils.FeatureDB(f'{db_paths["gffcompare"]}')
    reference_db = gffutils.FeatureDB(f'{db_paths["reference"]}')

    return gffcompare_db, reference_db
