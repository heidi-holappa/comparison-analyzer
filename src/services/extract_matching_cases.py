class MatchingCasesExtractor:

    def __init__(self, offset_results: dict, offset_range: tuple, reference_db):
        """
            Extract candidates matching the selected offset.

        Args:
            offset_results (dict): a dictionary containing the offset results for each transcript
            offset (int): the offset to match
            reference_db: reference to a gffutils database with the reference annotation
        """
        self.offset_results = offset_results
        self.offset_range = offset_range
        self.reference_db = reference_db

    def extract_candidates_matching_selected_offset(self) -> dict:
        extracted_candidates = {}
        for transcript_id, value in self.offset_results.items():
            strand, offsets = value["strand"], value["offsets"]
            reference_id = value["reference_id"]
            if strand == "-":
                offsets.reverse()
            for offset_exon_idx in range(1, len(offsets)-1):
                offset_exon_number = offset_exon_idx + 1
                if abs(offsets[offset_exon_idx][0]) >= self.offset_range[0] and \
                        abs(offsets[offset_exon_idx][0]) <= self.offset_range[1]:
                    entry_key = transcript_id + ".exon_" + \
                        str(offset_exon_number) + ".start" + \
                        '.offset_' + str(offsets[offset_exon_idx][0])
                    for exon in self.reference_db.children(reference_id, featuretype='exon'):
                        if int(exon['exon_number'][0]) == offset_exon_number:
                            extracted_candidates[entry_key] = {
                                "transcript_id": transcript_id,
                                "strand": strand,
                                "exon_number": offset_exon_number,
                                "location_type": "start",
                                "location": exon.start + offsets[offset_exon_idx][0],
                                "offset": offsets[offset_exon_idx][0]
                            }
                if abs(offsets[offset_exon_idx][1]) >= self.offset_range[0] and \
                        abs(offsets[offset_exon_idx][1]) <= self.offset_range[1]:
                    entry_key = transcript_id + ".exon_" + \
                        str(offset_exon_number) + ".end" + \
                        '.offset_' + str(offsets[offset_exon_idx][1])
                    for exon in self.reference_db.children(reference_id, featuretype='exon'):
                        if int(exon['exon_number'][0]) == offset_exon_number:
                            extracted_candidates[entry_key] = {
                                "transcript_id": transcript_id,
                                "strand": strand,
                                "exon_number": offset_exon_number,
                                "location_type": "end",
                                "location": exon.end + offsets[offset_exon_idx][1],
                                "offset": offsets[offset_exon_idx][1]
                            }
        return extracted_candidates
