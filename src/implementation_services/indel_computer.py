import pysam

from services.output_manager import default_output_manager as output_manager
from services.log_manager import default_log_manager as log_manager


def extract_location_from_cigar_string(cigar_tuples: list,
                                       reference_start: int,
                                       reference_end: int,
                                       location: int):
    relative_position = location - reference_start
    alignment_position = 0
    ref_position = 0

    for cigar_code in cigar_tuples:

        if cigar_code[0] in [0, 2, 3, 7, 8]:
            ref_position += cigar_code[1]
        if ref_position <= relative_position and not \
                reference_start + ref_position == reference_end:
            alignment_position += cigar_code[1]
        else:
            return alignment_position + (cigar_code[1] - (ref_position - relative_position))

    return -1


def count_indels_from_cigar_codes_in_given_window(cigar_tuples: list,
                                                  aligned_location: int,
                                                  direction: str,
                                                  indel_count: dict,
                                                  window_size: int):
    """
    Get cigar codes in a given window.

    Args:
        cigar_tuples (list): list of cigar tuples (cigar code, aligned position)
        aligned_location (int): aligned location
        loc_type (str): type of location (start or end)
    """
    result = {
        'deletions': 0,
        'insertions': 0
    }

    deletions = 0
    insertions = 0

    cigar_code_list = []
    location = 0

    debug_list = []

    if direction == "left":
        aligned_location = aligned_location - window_size + 1

    for cigar_code in cigar_tuples:
        if window_size == len(cigar_code_list):
            break
        if location + cigar_code[1] > aligned_location:
            overlap = location + \
                cigar_code[1] - (aligned_location + len(cigar_code_list))
            cigar_code_list.extend(
                [cigar_code[0] for _ in range(min(window_size -
                                                  len(cigar_code_list), overlap))])
        location += cigar_code[1]

    debug_list.append(cigar_code_list)
    # for cigar_code in cigar_code_list:
    #     if cigar_code == 2:
    #         deletions += 1
    #     if cigar_code == 1:
    #         insertions += 1

    for i in range(window_size):
        if i >= len(cigar_code_list):
            break
        if cigar_code_list[i] == 2:
            deletions += 1
            indel_count["del_pos_distr"][i] += 1
        if cigar_code_list[i] == 1:
            insertions += 1
            indel_count["ins_pos_distr"][i] += 1

    if deletions not in indel_count["deletions"]:
        indel_count["deletions"][deletions] = 0
    if insertions not in indel_count["insertions"]:
        indel_count["insertions"][insertions] = 0

    indel_count["deletions"][deletions] += 1
    indel_count["insertions"][insertions] += 1


def execute_indel_computation(
        bam_path: str,
        intron_site_dict: dict,
        reads_and_references: dict,
        window_size: int):

    output_manager.output_line({
        "line": "ISOQUANT-GTF: COUNTING INDELS",
        "is_title": True
    })
    output_manager.output_line({
        "line": "Iterating reads and counting indels. This may take a while.",
        "is_info": True
    })

    samfile = pysam.AlignmentFile(bam_path, "rb")

    count = 0
    errors = []
    set_of_processed_reads = set()
    prev_processed_reads_counter = 0
    for read in samfile.fetch():

        if read.is_supplementary:
            continue
        if read.query_name in reads_and_references:
            if read.query_name not in set_of_processed_reads:
                set_of_processed_reads.add(read.query_name)
            else:
                prev_processed_reads_counter += 1
                continue
            count += 1
            if count % 1000 == 0:
                output_manager.output_line({
                    "line": "Processed " + str(count) + " reads",
                    "end_line": "\r",
                    "is_info": True,
                    "save_to_log": False
                })
            for matching_case_key in reads_and_references[read.query_name]:
                location = intron_site_dict[matching_case_key]["location"]

                idx_corrected_location = location - 1

                if read.reference_start > idx_corrected_location \
                        or read.reference_end < idx_corrected_location:
                    errors.append(
                        f"Non-matching location: {read.query_name}, {matching_case_key}\t")
                    continue

                if not read.cigartuples:
                    continue

                if not read.reference_end:
                    errors.append(f"{read.query_name}\tno reference end\n")
                    continue

                for direction, collected_info in intron_site_dict[matching_case_key]["extracted_information"].items():
                    aligned_location = extract_location_from_cigar_string(
                        read.cigartuples,
                        read.reference_start,
                        read.reference_end,
                        idx_corrected_location
                    )

                    count_indels_from_cigar_codes_in_given_window(
                        read.cigartuples,
                        aligned_location,
                        direction,
                        collected_info,
                        window_size)

    output_manager.output_line({
        "line": "Processing BAM-file finished.",
        "is_info": True
    })

    if errors:
        log_manager.debug_logs["BAM-file-errors"] = errors
