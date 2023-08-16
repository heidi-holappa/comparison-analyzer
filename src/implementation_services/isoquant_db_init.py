import os
import gffutils
from services.output_manager import default_output_manager as output_manager


def init_isoquant_db(isoquant_gtf, force=False):
    gtf_paths = {
        'isoquant': isoquant_gtf
    }

    db_paths = {
        'isoquant': isoquant_gtf[:-4] + '-ca.db'
    }

    output_manager.output_line({
        "line": "INITIALIZING ISOQUANT GTF-DB",
        "is_title": True
    })
    output_manager.output_line({
        "line": "IsoQuant GTF-file: " + os.path.basename(isoquant_gtf),
        "is_info": True
    })

    for key, value in db_paths.items():
        db_exists = os.path.exists(f'{value}')

        if not force and db_exists:
            output_manager.output_line({
                "line": f"{key}: database file already exists. Using existing database. " +
                "Use -f to force overwriting existing db-files.",
                "is_info": True
            })

        else:
            output_manager.output_line({
                "line": f"{key}: database file does not exist. " +
                "Creating database. This might take a while.",
                "is_info": True
            })
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
            output_manager.output_line({
                "line": f"{key}: database created successfully!",
                "is_info": True
            })

    isoquant_db = gffutils.FeatureDB(f'{db_paths["isoquant"]}')

    return isoquant_db
