import os
import gffutils
import argparse

parser = argparse.ArgumentParser(
    description='compAna: a tool for comparing annotations',
    usage='python3 compAna.py -i <input.gtf> [-f] [-s]'
    )
parser.add_argument('-i', '--input', help='gtf-file to be imported into the database', required=True, metavar='')
parser.add_argument('-r', '--reference', help='reference gtf-file to be compared against', required=True, metavar='')
parser.add_argument('-f', '--force', help='force overwrite of existing database', action='store_true')
parser.add_argument('-s', '--stats', help='output statistics of class codes', action='store_true')
parser.add_argument('-c', '--class-code', help='specify gffcompare class code to analyze.', action='store_true')

arguments = parser.parse_args()

gtf_paths = {
    'gffcompare-gtf': arguments.input,
    'reference-gtf': arguments.reference
}

db_paths = {
    'gffcompare-db': arguments.input[:-4] + '-ca.db',
    'reference-db': arguments.reference[:-4] + '-ca.db'
}

print("=========================================")
print("compAna: a tool for comparing annotations")
print("=========================================")
print("===========DATABASE MANAGEMENT===========\n")

for key, value in db_paths.items():
    db_exists = os.path.exists(f'{value}')

    if not arguments.force and db_exists:
        print(f"{key}: using existing db file. Use -f to force overwrite existing db-files.")
        # gffutils_db = gffutils.FeatureDB(f'{db_path}/{db_name}')
    else:
        print(f'{key}: creating database... this might take a while.')
        gffutils.create_db(
            arguments.input, 
            dbfn=f'{value}', 
            force=True, 
            keep_order=True, 
            merge_strategy='merge', 
            sort_attribute_values=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True
            )
        print(f"{key}: database created successfully!")

gffcompare_db = gffutils.FeatureDB(f'{db_paths["gffcompare-db"]}')
reference_db = gffutils.FeatureDB(f'{db_paths["reference-db"]}')

print("\n=========================================\n")

if arguments.stats:
    print("==========CLASS CODE STATISTICS==========")
    class_codes = {}
    print('\nComputing statistics for class codes...\n')
    for feature in gffcompare_db.features_of_type('transcript'):
        if 'class_code' in feature.attributes:
            class_code = feature.attributes['class_code'][0]
            if not class_code in class_codes:
                class_codes[class_code] = 0
            class_codes[class_code] += 1

    print(f"{'class code':<18}| {'n':<10}")
    print('-' * 30)
    for key, value in sorted(class_codes.items(), key=lambda item: item[1], reverse=True):
        print(f"{key:<18}| {value:<10}")

print("\n=========================================")

if arguments.class_code:
    print("==========ANNOTATION COMPARISON==========")
    class_code = input("Enter class code to analyze: ")
    print(f"Analyzing class code: {class_code}")
    offset_analysis = {}
    for feature in gffcompare_db.features_of_type('transcript'):
        if 'class_code' in feature.attributes and class_code in feature.attributes['class_code']:
            for exon in gffcompare_db.children(feature, featuretype='exon', order_by='start'):
                if feature.id not in offset_analysis:
                    offset_analysis[feature.id] = []
                offset_analysis[feature.id].append((exon.start, exon.end))
    for key, value in offset_analysis.items():
        print(key, value)
                
    print("=========================================")




print("================= END ===================")
