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
            offset_between_aligned_and_reference_exon = calculate_total_offset(aligned_exons[e_index], reference_exons[r_index])
            if e_index < len(aligned_exons) - 1:
                offset_between_reference_exon_and_next_aligned_exon = calculate_total_offset(aligned_exons[e_index+1], reference_exons[r_index])
                if offset_between_aligned_and_reference_exon > offset_between_reference_exon_and_next_aligned_exon:
                    result = (float('inf'), float('inf'))
                    break
            if r_index < len(reference_exons) - 1:
                offset_between_aligned_exon_and_next_reference_exon = calculate_total_offset(aligned_exons[e_index], reference_exons[r_index+1])
                if offset_between_aligned_and_reference_exon > offset_between_aligned_exon_and_next_reference_exon:
                    if e_index < len(aligned_exons) - 1:
                        offset_between_next_reference_exon_and_next_aligned_exon = calculate_total_offset(aligned_exons[e_index+1], reference_exons[r_index+1])
                        if not offset_between_next_reference_exon_and_next_aligned_exon < offset_between_aligned_exon_and_next_reference_exon:
                            offset_list.append((float('-inf'), float('-inf')))
                            r_start_index = r_index + 1
                            continue
                    else:
                        offset_list.append((float('-inf'), float('-inf')))
                        r_start_index = r_index + 1
                        continue
            result = (aligned_exons[e_index][0] - reference_exons[r_index][0], aligned_exons[e_index][1] - reference_exons[r_index][1])
            r_start_index = r_index + 1
            break
        offset_list.append(result)
    if r_start_index < len(reference_exons):
        for r_index in range(r_start_index, len(reference_exons)):
            offset_list.append((float('-inf'), float('-inf')))
    return offset_list