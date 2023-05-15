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
parser.add_argument('-c', '--class-code', nargs='+', help='specify gffcompare class code to analyze.')

arguments = parser.parse_args()

gtf_paths = {
    'gffcompare': arguments.input,
    'reference': arguments.reference
}

db_paths = {
    'gffcompare': arguments.input[:-4] + '-ca.db',
    'reference': arguments.reference[:-4] + '-ca.db'
}

print("=========================================")
print("compAna: a tool for comparing annotations")
print("=========================================")
print("===========DATABASE MANAGEMENT===========\n")

print("============ FILE INFORMATION ===========\n")
print(f"Gffcompare GTF-file: {os.path.basename(arguments.input)}")
print(f"Reference GTF-file: {os.path.basename(arguments.reference)}\n")


for key, value in db_paths.items():
    db_exists = os.path.exists(f'{value}')

    if not arguments.force and db_exists:
        print(f"{key}: using existing db file. Use -f to force overwrite existing db-files.")
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

if arguments.stats:
    print("==========CLASS CODE STATISTICS==========")
    class_codes = {}
    print('\nComputing statistics for class codes...\n')
    for transcript in gffcompare_db.features_of_type('transcript'):
        if 'class_code' in transcript.attributes:
            class_code = transcript.attributes['class_code'][0]
            if not class_code in class_codes:
                class_codes[class_code] = 0
            class_codes[class_code] += 1

    print(f"{'class code':<18}| {'n':<10}")
    print('-' * 30)
    for key, value in sorted(class_codes.items(), key=lambda item: item[1], reverse=True):
        print(f"{key:<18}| {value:<10}")

print("\n=========================================")

if arguments.class_code:
    
    for class_code in arguments.class_code:
        print("==========ANNOTATION COMPARISON==========")
        print(f"Analyzing class code: {class_code}")
        offset_analysis = {}
        for transcript in gffcompare_db.features_of_type('transcript'):
            if 'class_code' in transcript.attributes and class_code in transcript.attributes['class_code']:
                ref_transcript_id = transcript['cmp_ref'][0]
                for exon in gffcompare_db.children(transcript, featuretype='exon', order_by='start'):
                    for reference_exon in reference_db.children(ref_transcript_id, featuretype='exon', order_by='start'):
                        if transcript['cmp_ref'] != reference_exon['transcript_id']:
                            continue
                        dict_key = (transcript.id, reference_exon['transcript_id'][0])
                        if dict_key not in offset_analysis:
                            offset_analysis[dict_key] = {}
                        offset = (exon.start - reference_exon.start, exon.end - reference_exon.end)
                        exon_number = int(exon['exon_number'][0])
                        if exon_number not in offset_analysis[dict_key]:
                            offset_analysis[dict_key][exon_number] = offset
                        elif abs(offset[1]) + abs(offset[0]) < abs(offset_analysis[dict_key][exon_number][0]) + (offset_analysis[dict_key][exon_number][1]):
                            offset_analysis[dict_key][exon_number] = offset
        for key, value in offset_analysis.items():
            print(key, value)
                    
        print("=========================================\n")




print("================= END ===================")
