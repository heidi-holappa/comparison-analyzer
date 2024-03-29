from services.output_manager import default_output_manager as output_manager
from services.log_manager import default_log_manager as log_manager


def calculate_total_offset(exon_1, exon_2):
    start_offset = exon_1[0] - exon_2[0]
    end_offset = exon_1[1] - exon_2[1]
    total_offset = abs(start_offset) + abs(end_offset)
    return total_offset


def compute_offsets(aligned_exons, reference_exons):
    r_start_index = 0
    offset_list = []
    for e_index in range(0, len(aligned_exons)):
        result = (float('inf'), float('inf'))
        for r_index in range(r_start_index, len(reference_exons)):
            offset_aligned_and_reference_exon = calculate_total_offset(
                aligned_exons[e_index], reference_exons[r_index])
            if e_index < len(aligned_exons) - 1:
                offset_reference_exon_and_next_aligned_exon = calculate_total_offset(
                    aligned_exons[e_index+1], reference_exons[r_index])
                if offset_aligned_and_reference_exon > offset_reference_exon_and_next_aligned_exon:
                    result = (float('inf'), float('inf'))
                    break
            if r_index < len(reference_exons) - 1:
                offset_aligned_exon_and_next_reference_exon = calculate_total_offset(
                    aligned_exons[e_index], reference_exons[r_index+1])
                if offset_aligned_and_reference_exon > offset_aligned_exon_and_next_reference_exon:
                    if e_index < len(aligned_exons) - 1:
                        offset_next_reference_exon_and_next_aligned_exon = calculate_total_offset(
                            aligned_exons[e_index+1], reference_exons[r_index+1])
                        if not offset_next_reference_exon_and_next_aligned_exon \
                                < offset_aligned_exon_and_next_reference_exon:
                            offset_list.append((float('-inf'), float('-inf')))
                            r_start_index = r_index + 1
                            continue
                    else:
                        offset_list.append((float('-inf'), float('-inf')))
                        r_start_index = r_index + 1
                        continue
            result = (aligned_exons[e_index][0] - reference_exons[r_index]
                      [0], aligned_exons[e_index][1] - reference_exons[r_index][1])
            r_start_index = r_index + 1
            break
        offset_list.append(result)
    if r_start_index < len(reference_exons):
        for r_index in range(r_start_index, len(reference_exons)):
            offset_list.append((float('-inf'), float('-inf')))
    return offset_list


def fetch_exons(transcript, class_code, gffcompare_db, reference_db):
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
    if not 'class_code' in transcript.attributes or \
            not class_code in transcript.attributes['class_code']:
        return aligned_exons, reference_exons
    for exon in gffcompare_db.children(transcript, featuretype='exon', order_by='start'):
        aligned_exons.append((exon.start, exon.end))
    ref_transcript_id = transcript['cmp_ref'][0]
    for exon in reference_db.children(ref_transcript_id, featuretype='exon', order_by='start'):
        if transcript['cmp_ref'] != exon['transcript_id']:
            continue
        reference_exons.append((exon.start, exon.end))
    return aligned_exons, reference_exons


def execute_offset_computation(
        class_code: str,
        gffcompare_db,
        reference_db,
        extended_debug: bool) -> dict:
    output_manager.output_line({
        "line": "COMPUTE OFFSETS",
        "is_title": True
    })
    output_manager.output_line({
        "line": "Computing offsets for Gffcompare class codes: " + class_code,
        "is_info": True
    })

    offset_results_as_dict = {}
    for class_code in class_code.strip().split(" "):

        for transcript_item in gffcompare_db.features_of_type('transcript'):
            aligned_exons, reference_exons = fetch_exons(
                transcript_item,
                class_code,
                gffcompare_db,
                reference_db
            )
            if aligned_exons:
                offsets = compute_offsets(aligned_exons, reference_exons)
                offset_results_as_dict[transcript_item.id] = {
                    "offsets": offsets,
                    "reference_id": transcript_item['cmp_ref'][0],
                    "strand": transcript_item.strand,
                    "class_code": class_code
                }

    if extended_debug:
        log_manager.offset_results = offset_results_as_dict
        output_manager.output_line({
            "line": "Finished.",
            "is_info": True
        })
    else:
        output_manager.output_line({
            "line": "Offset computation finished. " +
            "Hint: to output results into a file, enable extended debug.",
            "is_info": True
        })
    return offset_results_as_dict
