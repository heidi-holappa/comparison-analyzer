import os
import json
import gffutils
import argparse
from services.offset_computation import compute_offsets
from services.fasta_extractor import FastaExtractor


parser = argparse.ArgumentParser(
    description='compAna: a tool for comparing annotations',
    usage='python3 compAna.py -i <input.gtf> [-f] [-s]'
    )
parser.add_argument('-g', '--gffcompare_gtf', help='gtf-file to be imported into the database', required=False, metavar='')
parser.add_argument('-r', '--reference_gtf', help='reference gtf-file to be compared against', required=False, metavar='')
parser.add_argument('-a', '--reference_fasta', help='reference fasta-file to be compared against', required=False, metavar='')
parser.add_argument('-o', '--offset', help='offset from which canonical splice point is to be searched', required=False, metavar='')
parser.add_argument('-f', '--force', help='force overwrite of existing database', action='store_true')
parser.add_argument('-s', '--stats', help='output statistics of class codes', action='store_true')
parser.add_argument('-c', '--class-code', nargs='+', help='specify gffcompare class code to analyze.')
parser.add_argument('-j', '--json', help='input arguments from json file', metavar='')

arguments = parser.parse_args()
argparse_dict = vars(arguments)


if arguments.json:
    with open(arguments.json) as json_file:
        json_dict = json.load(json_file)
        argparse_dict.update(json_dict)
        print(argparse_dict)


    


gtf_paths = {
    'gffcompare': arguments.gffcompare_gtf,
    'reference': arguments.reference_gtf
}

db_paths = {
    'gffcompare': arguments.gffcompare_gtf[:-4] + '-ca.db',
    'reference': arguments.reference_gtf[:-4] + '-ca.db'
}

print("=========================================")
print("compAna: a tool for comparing annotations")
print("=========================================")
print("===========DATABASE MANAGEMENT===========\n")

print("============ FILE INFORMATION ===========\n")
print(f"Gffcompare GTF-file: {os.path.basename(arguments.gffcompare_gtf)}")
print(f"Reference GTF-file: {os.path.basename(arguments.reference_gtf)}\n")


for key, value in db_paths.items():
    """
        Check if database files exists. If they do, use existing files. If not, create db-file(s).
    """
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
    """
        Compute simple n-count statistics for class codes.
    """ 
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
    """
        Fetch exons from gffcompare and reference databases.

    Args:
        transcript (gffutils.feature.Feature): a feature of type 'transcript'
        class_code (str): a character representing the class code of the transcript

    Returns:
        _type_: _description_
    """
    aligned_exons = []
    reference_exons = []
    if not 'class_code' in transcript.attributes or not class_code in transcript.attributes['class_code']:
        return aligned_exons, reference_exons
    for exon in gffcompare_db.children(transcript, featuretype='exon', order_by='start'):
        aligned_exons.append((exon.start, exon.end))
    ref_transcript_id = transcript['cmp_ref'][0]
    for exon in reference_db.children(ref_transcript_id, featuretype='exon', order_by='start'):
        if transcript['cmp_ref'] != exon['transcript_id']:
            continue
        reference_exons.append((exon.start, exon.end))
    return aligned_exons, reference_exons



if arguments.class_code:
    """
        Compute offsets for specified class codes.
    """
    offset_results = {}
    for class_code in arguments.class_code:
        print("==========ANNOTATION COMPARISON==========")
        print(f"Analyzing class code: {class_code}")
        class_code_results = {}
        for transcript in gffcompare_db.features_of_type('transcript'):
            aligned_exons, reference_exons = fetch_exons(transcript, class_code)
            if aligned_exons:
                offsets = compute_offsets(aligned_exons, reference_exons)
                dict_key = (transcript.id, transcript['cmp_ref'][0], transcript.strand)
                class_code_results[dict_key] = offsets
                offset_results[dict_key] = offsets
        for key, value in class_code_results.items():
            print(f"{key}: {value}")
                    
        print("=========================================\n")

def extract_candidates_matching_selected_offset(offset_results: dict, offset: int):
    """
        Extract candidates matching the selected offset.

    Args:
        offset_results (dict): a dictionary containing the offset results for each transcript
        offset (int): the offset to match
    """
    extracted_candidates = {}
    for key, value in offset_results.items():
        for i in range(1, len(value)-1):
            if abs(value[i][0]) == offset:
                for exon in reference_db.children(key[1], featuretype='exon', order_by='start'):
                    if int(exon['exon_number'][0]) == i + 1:
                        extracted_candidates[(key[0], key[1], key[2], i + 1, 'start')] = exon.start
                        break
            elif abs(value[i][1]) == offset: 
                for exon in reference_db.children(key[1], featuretype='exon', order_by='start'):
                    if int(exon['exon_number'][0]) == i + 1:  
                        extracted_candidates[(key[0], key[1], key[2], i + 1, 'end')] = exon.end
                        break
    return extracted_candidates
                


print("==========CHARACTERS AT OFFSET===========\n")
if arguments.reference_fasta:
    print("Fetching reference fasta file...", end=' ')
    try:
        fasta_extractor = FastaExtractor(arguments.reference_fasta)
        print("success!\n")
    except:
        print("fasta-file not found. Please check path and try again.")
    if not arguments.offset:
        print("No offset value given. Nothing to do here.")
    else:
        print(f"Extracting candidates matching offset {arguments.offset}...")
        extracted_candidates = extract_candidates_matching_selected_offset(offset_results, arguments.offset)
        # for key, value in extracted_candidates.items():
        #     print(f"{key}: {value}")    
        results = {}
        for key, value in extracted_candidates.items():
            chromosome = key[0].split('.')[1]
            if key[4] == 'end':
                coordinates = (chromosome, value, value + 2)
            else:
                coordinates = (chromosome, value - 3, value - 1)
            chars = fasta_extractor.extract_characters_at_given_coordinates(coordinates)
            results[key] = chars
        json_overview = {
            "strand": {
                "+": {
                    "start": {},
                    "end": {}
                },
                "-": {
                    "start": {},
                    "end": {}
                }
            },
            
        }
        for key, value in results.items():
            if str(value) not in json_overview['strand'][key[2]][key[4]]:
                json_overview['strand'][key[2]][key[4]][str(value)] = 0
            json_overview['strand'][key[2]][key[4]][str(value)] += 1
        
        for strand in json_overview['strand']:
            for position in json_overview['strand'][strand]:
                json_overview['strand'][strand][position] = dict(sorted(json_overview['strand'][strand][position].items(), key=lambda item: item[1], reverse=True))


        with open("overview.md", "w") as file:
            file.write("# Overview\n")
            file.write("## Offset characters\n")
            file.write("This section contains the characters in the reference FASTA-file at the offset specified by the user.  \n\n")
            file.write("**Arguments provided by the user:**\n")
            file.write("```\n")
            file.write("gffcompare GTF-file:\n" + arguments.gffcompare_gtf + "\n\n")
            file.write("Reference GTF-file:\n" + arguments.reference_gtf + "\n\n")
            file.write("Reference FASTA-file:\n" + arguments.reference_fasta + "\n\n")
            file.write("Specified offset: " + str(arguments.offset) + "\n")
            file.write("Class codes: " + str(arguments.class_code) + "\n")
            file.write("```\n")
            file.write("**Results in JSON-format:**  \n")
            file.write("```json\n")
            file.write(json.dumps(json_overview, indent=4))
            file.write("\n```\n")
            

        
    




print("================= END ===================")
