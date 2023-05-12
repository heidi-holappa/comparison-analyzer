import os
import gffutils
import argparse

parser = argparse.ArgumentParser(
    description='compAna: a tool for comparing annotations',
    usage='python3 compAna.py -i <input.gtf> [-f] [-s]'
    )
parser.add_argument('-i', '--input', help='gtf-file to be imported into the database', required=True, metavar='')
parser.add_argument('-r', '--reference', help='reference gtf-file to be compared against', metavar='')
parser.add_argument('-f', '--force', help='force overwrite of existing database', action='store_true')
parser.add_argument('-s', '--stats', help='output statistics of class codes', action='store_true')

arguments = parser.parse_args()


db_path = os.path.dirname(arguments.input)
db_name = os.path.basename(arguments.input)[:-4] + '-ca.db'
db_exists = os.path.exists(f'{db_path}/{db_name}')

if not arguments.force and db_exists:
    print("Using existing db file. Use -f to force overwrite.")
    gffutils_db = gffutils.FeatureDB(f'{db_path}/{db_name}')
else:
    print('Creating database... this might take a while.')
    gffutils_db = gffutils.create_db(
        arguments.input, 
        dbfn=f'{db_path}/{db_name}', 
        force=True, 
        keep_order=True, 
        merge_strategy='merge', 
        sort_attribute_values=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True
        )
    print("Database created successfully!")

if arguments.stats:
    class_codes = {}
    print('Computing statistics for class codes...')
    for feature in gffutils_db.all_features():
        if 'class_code' in feature.attributes:
            class_code = feature.attributes['class_code'][0]
            if not class_code in class_codes:
                class_codes[class_code] = 0
            class_codes[class_code] += 1

    print(f"{'class code':<18}| {'n':<10}")
    print('-' * 30)
    for key, value in sorted(class_codes.items(), key=lambda item: item[1], reverse=True):
        print(f"{key:<18}| {value:<10}")



# Sample file path
# /home/holaphei/koulutyot/lv-2022-2023/BI-summer-trainee/isoquant-summer-trainee/ca-task-2/comparison-analyzer/sample-data/nanopore-compare-2.annotated.gtf