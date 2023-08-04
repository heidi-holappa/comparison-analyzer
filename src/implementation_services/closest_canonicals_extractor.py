from pyfaidx import Fasta


def initialize_fasta(fasta_path):
    return Fasta(fasta_path)


def extract_characters_at_given_coordinates(coordinates: tuple, index_correction: int, fasta: Fasta):
    chromosome, start, end = coordinates
    return fasta[chromosome][start + index_correction:end + index_correction]


def find_closest_canonicals(nucleotides: str, dict_key: str, canonicals: list, intron_site_dict: dict):
    nucleotides_middle = int(len(nucleotides) / 2)
    closest_canonicals = {}
    aligned_splice_site_nucleotides = nucleotides[nucleotides_middle:nucleotides_middle + 2]
    for i in range(1, nucleotides_middle):
        left_position = nucleotides[nucleotides_middle -
                                    i:nucleotides_middle - i + 2]
        right_position = nucleotides[nucleotides_middle +
                                     i:nucleotides_middle + i + 2]
        if left_position in canonicals and 'left' not in closest_canonicals:
            closest_canonicals['left'] = (
                left_position, aligned_splice_site_nucleotides, i)
        if right_position in canonicals and 'right' not in closest_canonicals:
            closest_canonicals['right'] = (
                right_position, aligned_splice_site_nucleotides, i)
    for item in ['left', 'right']:
        if item not in closest_canonicals:
            closest_canonicals[item] = (
                aligned_splice_site_nucleotides, aligned_splice_site_nucleotides, 0)
    intron_site_dict[dict_key]["extracted_information"]["left"]['closest_canonical'] = closest_canonicals['left']
    intron_site_dict[dict_key]["extracted_information"]["right"]['closest_canonical'] = closest_canonicals['right']


def iterate_intron_sites(intron_site_dict: dict, window_size: int, index_correction: int, fasta: Fasta):
    for key, value in intron_site_dict.items():
        # TODO: check if this is correct and get rid of magic numbers
        if value["location_type"] == "start":
            splice_cite_location = value["location"] - 2
        else:
            splice_cite_location = value["location"] + 1
        coordinates = (value["seq_id"],
                       splice_cite_location - window_size,
                       splice_cite_location + window_size)
        nucleotides = extract_characters_at_given_coordinates(
            coordinates, index_correction, fasta)

        possible_canonicals = {
            '+': {
                'start': ['AG', 'AC'],
                'end': ['GT', 'GC', 'AT']
            },
            '-': {
                'start': ['AC', 'GC', 'AC'],
                'end': ['CT', 'GT']
            }
        }
        strand = value["strand"]
        location_type = value["location_type"]
        canonicals = possible_canonicals[strand][location_type]

        find_closest_canonicals(str(nucleotides), key,
                                canonicals, intron_site_dict)


def execute_closest_canonicals_extraction(intron_site_dict: dict, window_size: int, fasta_path: str):
    index_correction = -1
    fasta = initialize_fasta(fasta_path)
    iterate_intron_sites(intron_site_dict, window_size,
                         index_correction, fasta)
