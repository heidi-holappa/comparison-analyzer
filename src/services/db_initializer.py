import os
import gffutils


def init_databases(parser):
    gtf_paths = {
        'gffcompare': parser.gffcompare_gtf,
        'reference': parser.reference_gtf
    }

    db_paths = {
        'gffcompare': parser.gffcompare_gtf[:-4] + '-ca.db',
        'reference': parser.reference_gtf[:-4] + '-ca.db'
    }

    print("============ FILE INFORMATION ===========\n")
    print(f"Gffcompare GTF-file: {os.path.basename(parser.gffcompare_gtf)}")
    print(f"Reference GTF-file: {os.path.basename(parser.reference_gtf)}\n")

    for key, value in db_paths.items():
        db_exists = os.path.exists(f'{value}')

        if not parser.force and db_exists:
            print(
                f"{key}: using existing db file. Use -f to force overwrite existing db-files.")
            # gffutils_db = gffutils.FeatureDB(f'{db_path}/{db_name}')
        else:
            print(f'{key}: creating database... this might take a while.')
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
            print(f"{key}: database created successfully!")

    gffcompare_db = gffutils.FeatureDB(f'{db_paths["gffcompare"]}')
    reference_db = gffutils.FeatureDB(f'{db_paths["reference"]}')

    print("\n=========================================\n")
    return gffcompare_db, reference_db
