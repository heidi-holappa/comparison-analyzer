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

def fetch_exons(transcript, class_code):
    analyzed_exons = []
    reference_exons = []
    if not 'class_code' in transcript.attributes or not class_code in transcript.attributes['class_code']:
        return analyzed_exons, reference_exons
    for exon in gffcompare_db.children(transcript, featuretype='exon', order_by='start'):
        analyzed_exons.append((exon.start, exon.end))
    ref_transcript_id = transcript['cmp_ref'][0]
    for exon in reference_db.children(ref_transcript_id, featuretype='exon', order_by='start'):
        if transcript['cmp_ref'] != exon['transcript_id']:
            continue
        reference_exons.append((exon.start, exon.end))
    return analyzed_exons, reference_exons

def calculate_total_offset(exon_1, exon_2):
    start_offset = exon_1[0] - exon_2[0]
    end_offset = exon_1[1] - exon_2[1]
    total_offset = abs(start_offset) + abs(end_offset)
    return total_offset

def compute_offset(analyzed_exons, reference_exons):
    r_start_index = 0
    offset_list = []
    for e_index in range(0, len(analyzed_exons)):
        result = (float('inf'), float('inf'))
        for r_index in range(r_start_index, len(reference_exons)):
            total_offset_current_analysis_exon_to_current_reference_exon = calculate_total_offset(analyzed_exons[e_index], reference_exons[r_index])
            if e_index < len(analyzed_exons) - 1:
                total_offset_current_reference_exon_to_next_analysis_exon = calculate_total_offset(analyzed_exons[e_index+1], reference_exons[r_index])
                if r_index < len(reference_exons) - 1 and total_offset_current_reference_exon_to_next_analysis_exon < total_offset_current_analysis_exon_to_current_reference_exon:
                    total_offset_from_next_analysis_exon_to_next_reference_exon = calculate_total_offset(analyzed_exons[e_index+1], reference_exons[r_index+1])
                    if total_offset_current_reference_exon_to_next_analysis_exon < total_offset_from_next_analysis_exon_to_next_reference_exon:
                        r_start_index = r_index
                        break
            if total_offset_current_analysis_exon_to_current_reference_exon < abs(result[0]) + abs(result[1]):
                if result != (float('inf'), float('inf')):
                    offset_list.append( (float('-inf'), float('-inf')))
                result = (analyzed_exons[e_index][0] - reference_exons[r_index][0], analyzed_exons[e_index][1] - reference_exons[r_index][1])
                r_start_index = r_index + 1
            else:
                break
        offset_list.append(result)
    return offset_list

if arguments.class_code:
        
    for class_code in arguments.class_code:
        print("==========ANNOTATION COMPARISON==========")
        print(f"Analyzing class code: {class_code}")
        offset_results = {}
        for transcript in gffcompare_db.features_of_type('transcript'):
            analyzed_exons, reference_exons = fetch_exons(transcript, class_code)
            if analyzed_exons:
                offsets = compute_offset(analyzed_exons, reference_exons)
                dict_key = (transcript.id, transcript['cmp_ref'][0], transcript.strand)
                offset_results[dict_key] = offsets
        for key, value in offset_results.items():
            print(f"{key}: {value}")
                    
        print("=========================================\n")



print("================= END ===================")
