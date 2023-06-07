import json

from services.fasta_extractor import FastaExtractor
from services.bam_manager import BamManager

from services.offset_computation import execute_offset_computation
from services.argument_parser import init_argparser
from services.db_initializer import init_databases
from services.class_code_stats import compute_class_code_stats

print("=========================================")
print("compAna: a tool for comparing annotations")
print("=========================================")
print("===========DATABASE MANAGEMENT===========\n")

parser = init_argparser()
gffcompare_db, reference_db = init_databases(parser)

if parser.stats:
    compute_class_code_stats(gffcompare_db)

if parser.class_code:
    execute_offset_computation(parser, gffcompare_db, reference_db)


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
                        extracted_candidates[(
                            key[0], key[1], key[2], i + 1, 'start')] = exon.start
                        break
            elif abs(value[i][1]) == offset:
                for exon in reference_db.children(key[1], featuretype='exon', order_by='start'):
                    if int(exon['exon_number'][0]) == i + 1:
                        extracted_candidates[(
                            key[0], key[1], key[2], i + 1, 'end')] = exon.end
                        break
    return extracted_candidates


print("==========CHARACTERS AT OFFSET===========\n")
if parser.reference_fasta:
    print("Fetching reference fasta file...", end=' ')
    try:
        fasta_extractor = FastaExtractor(parser.reference_fasta)
        print("success!\n")
    except:
        print("fasta-file not found. Please check path and try again.")
    if not parser.offset:
        print("No offset value given. Nothing to do here.")
    else:
        print(f"Extracting candidates matching offset {parser.offset}...")
        matching_cases_dict = extract_candidates_matching_selected_offset(
            offset_results, parser.offset)
        # for key, value in extracted_candidates.items():
        #     print(f"{key}: {value}")
        results = {}
        for key, value in matching_cases_dict.items():
            chromosome = key[0].split('.')[1]
            if key[4] == 'end':
                coordinates = (chromosome, value, value + 2)
            else:
                coordinates = (chromosome, value - 3, value - 1)
            chars = fasta_extractor.extract_characters_at_given_coordinates(
                coordinates)
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
                json_overview['strand'][strand][position] = dict(
                    sorted(
                        json_overview['strand'][strand][position].items(),
                        key=lambda item: item[1],
                        reverse=True
                    )
                )

        with open("overview.md", "w") as file:
            file.write("# Overview\n")
            file.write("## Offset characters\n")
            file.write(
                "This section contains the characters in the reference ")
            file.write(
                "FASTA-file at the offset specified by the user.  \n\n")
            file.write("**Arguments provided by the user:**\n")
            file.write("```\n")
            file.write("gffcompare GTF-file:\n" +
                       parser.gffcompare_gtf + "\n\n")
            file.write("Reference GTF-file:\n" +
                       parser.reference_gtf + "\n\n")
            file.write("Reference FASTA-file:\n" +
                       parser.reference_fasta + "\n\n")
            file.write("Specified offset: " + str(parser.offset) + "\n")
            file.write("Class codes: " + str(parser.class_code) + "\n")
            file.write("```\n")
            file.write("**Results in JSON-format:**  \n")
            file.write("```json\n")
            file.write(json.dumps(json_overview, indent=4))
            file.write("\n```\n")

if parser.reads_tsv and parser.reads_bam:
    print("Fetching reads from BAM-file. This might take some time.")
    transcripts = set()
    for row in matching_cases_dict:
        transcripts.add(row[0])
    bam_manager = BamManager(
        parser.reads_bam, parser.reads_tsv, transcripts)
    bam_manager.execute()


print("================= END ===================")
