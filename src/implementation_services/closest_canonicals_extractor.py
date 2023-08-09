from pyfaidx import Fasta


def initialize_fasta(fasta_path):
    return Fasta(fasta_path)


def extract_characters_at_given_coordinates(coordinates: tuple, index_correction: int, fasta: Fasta):
    chromosome, start, end = coordinates
    try:
        return fasta[chromosome][start + index_correction:end + index_correction]
    except KeyError:
        return "X" * (end - start)


def count_most_common_indel_case(indel_dict: dict):
    most_common_case = [
        k for k, v in indel_dict.items() if v == max(indel_dict.values())]
    if len(most_common_case) > 1 or len(most_common_case) == 0:
        return -1
    return most_common_case[0]


def extract_nucleotides_from_most_common_del_location(dict_entry: dict, fasta: Fasta):
    deletions_right = dict_entry["extracted_information"]["right"]["deletions"]
    deletions_left = dict_entry["extracted_information"]["left"]["deletions"]
    most_common_right = count_most_common_indel_case(deletions_right)
    most_common_left = count_most_common_indel_case(deletions_left)
    if most_common_right != -1:
        chr = dict_entry["seq_id"]
        location = dict_entry["location"]
        coordinates = (chr, location + most_common_right,
                       location + most_common_right + 2)
        dict_entry["extracted_information"]["right"]["most_common_case_nucleotides"] = str(extract_characters_at_given_coordinates(
            coordinates, -1, fasta))
    if most_common_left != -1:
        chr = dict_entry["seq_id"]
        location = dict_entry["location"]
        coordinates = (chr, location - most_common_left - 2,
                       location - most_common_left)
        dict_entry["extracted_information"]["left"]["most_common_case_nucleotides"] = str(extract_characters_at_given_coordinates(
            coordinates, -1, fasta))


def find_most_common_del_nucleotide_pair(nucleotides: str, dict_key: str, canonicals: list, intron_site_dict: dict):
    nucleotides_middle = int(len(nucleotides) / 2)
    # closest_canonicals = {}
    # aligned_splice_site_nucleotides = nucleotides[nucleotides_middle:nucleotides_middle + 2]
    deletions_right = intron_site_dict[dict_key]["extracted_information"]["right"]["deletions"]
    right_offset = count_most_common_indel_case(deletions_right)
    deletions_left = intron_site_dict[dict_key]["extracted_information"]["left"]["deletions"]
    left_offset = count_most_common_indel_case(deletions_left)
    if right_offset != -1:
        nucleotide_pair_right = nucleotides[nucleotides_middle +
                                            right_offset:nucleotides_middle + right_offset + 2]
    else:
        nucleotide_pair_right = "XX"
    if left_offset != -1:
        nucleotide_pair_left = nucleotides[nucleotides_middle -
                                           left_offset:nucleotides_middle - left_offset + 2]
    else:
        nucleotide_pair_left = "XX"

    intron_site_dict[dict_key]["extracted_information"]["left"]['most_common_del_pair'] = nucleotide_pair_left
    intron_site_dict[dict_key]["extracted_information"]["right"]['most_common_del_pair'] = nucleotide_pair_right
    dict_entry["extracted_information"]["right"]["most_common_del"] = deletions_right
    dict_entry["extracted_information"]["left"]["most_common_del"] = deletions_left


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
    intron_site_dict[dict_key]["nucleotides"] = nucleotides


def iterate_intron_sites(intron_site_dict: dict, window_size: int, index_correction: int, fasta: Fasta):
    for key, value in intron_site_dict.items():
        if value['strand'] not in ['+', '-']:
            intron_site_dict[key]["extracted_information"]["left"]['closest_canonical'] = (
                'XX', 'XX', 0)
            intron_site_dict[key]["extracted_information"]["right"]['closest_canonical'] = (
                'XX', 'XX', 0)
            continue

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

        # extract_nucleotides_from_most_common_del_location(value, fasta)

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

        find_most_common_del_nucleotide_pair(str(nucleotides), key,
                                             canonicals, intron_site_dict)


def execute_closest_canonicals_extraction(intron_site_dict: dict, window_size: int, fasta_path: str):
    index_correction = -1
    fasta = initialize_fasta(fasta_path)
    iterate_intron_sites(intron_site_dict, window_size,
                         index_correction, fasta)
