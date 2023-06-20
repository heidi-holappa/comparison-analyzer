class MatchingCasesExtractor:

    def __init__(self, offset_results: dict, offset: int, reference_db):
        """
            Extract candidates matching the selected offset.

        Args:
            offset_results (dict): a dictionary containing the offset results for each transcript
            offset (int): the offset to match
            reference_db: reference to a gffutils database with the reference annotation
        """
        self.offset_results = offset_results
        self.offset = offset
        self.reference_db = reference_db

    def extract_candidates_matching_selected_offset(self) -> dict:
        extracted_candidates = {}
        for key, value in self.offset_results.items():
            strand, offsets = value["strand"], value["offsets"]
            reference_id, transcript_id = value["reference_id"], key
            if strand == "-":
                offsets.reverse()
            for offset_exon_idx in range(1, len(offsets)-1):
                offset_exon_number = offset_exon_idx + 1
                if abs(offsets[offset_exon_idx][0]) == self.offset:
                    entry_key = transcript_id + ".exon_" + \
                        str(offset_exon_number) + ".start"
                    for exon in self.reference_db.children(reference_id, featuretype='exon'):
                        if int(exon['exon_number'][0]) == offset_exon_number:
                            extracted_candidates[entry_key] = {
                                "transcript_id": transcript_id,
                                "strand": strand,
                                "exon_number": offset_exon_number,
                                "location_type": "start",
                                "location": exon.start + offsets[offset_exon_idx][0]
                            }
                if abs(offsets[offset_exon_idx][1]) == self.offset:
                    entry_key = transcript_id + ".exon_" + \
                        str(offset_exon_number) + ".end"
                    for exon in self.reference_db.children(reference_id, featuretype='exon'):
                        if int(exon['exon_number'][0]) == offset_exon_number:
                            extracted_candidates[entry_key] = {
                                "transcript_id": transcript_id,
                                "strand": strand,
                                "exon_number": offset_exon_number,
                                "location_type": "end",
                                "location": exon.end + offsets[offset_exon_idx][1]
                            }
        return extracted_candidates
