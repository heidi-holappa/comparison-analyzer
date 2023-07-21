from services.output_manager import default_output_manager as output_manager


class MatchingCasesExtractor:

    def __init__(self, offset_results: dict, offset_range: tuple, reference_db):
        """
            Extract candidates matching the selected offset. Currently first and last 
            exon are not considered. These are omitted with the magic numbers in the 
            second for loop.

        Args:
            offset_results (dict): a dictionary containing the offset results for each transcript
            offset (int): the offset to match
            reference_db: reference to a gffutils database with the reference annotation
        """
        self.offset_results = offset_results
        self.offset_range = offset_range
        self.reference_db = reference_db

    def extract_candidates_matching_selected_offset(self) -> dict:
        output_manager.output_line({
            "line": "EXTRACTING CASES",
            "is_title": True
        })
        output_manager.output_line({
            "line": f"Extracting candidates in the given offset range {self.offset_range}",
            "is_info": True
        })
        extracted_candidates = {}
        for transcript_id, value in self.offset_results.items():
            strand, offsets, reference_id = value["strand"], value["offsets"], value["reference_id"]
            for offset_exon_idx in range(1, len(offsets)-1):
                start_is_in_range = bool(
                    abs(offsets[offset_exon_idx][0]) >= self.offset_range[0] and
                    abs(offsets[offset_exon_idx][0]) <= self.offset_range[1])
                if start_is_in_range:
                    entry_key = transcript_id + ".exon_" + \
                        str(offset_exon_idx) + ".start" + \
                        '.offset_' + str(offsets[offset_exon_idx][0])
                    for exon_number, exon in enumerate(self.reference_db.children(
                            reference_id, featuretype='exon', order_by='start')):
                        new_entry_key = (transcript_id, exon.start)
                        if exon_number == offset_exon_idx:
                            extracted_candidates[new_entry_key] = {
                                "transcript_id": transcript_id,
                                "strand": strand,
                                "exon_number": offset_exon_idx,
                                "location_type": "start",
                                "location": exon.start + offsets[offset_exon_idx][0],
                                "offset": offsets[offset_exon_idx][0],
                                "seq_id": exon.seqid
                            }
                end_is_in_range = bool(abs(offsets[offset_exon_idx][1]) >= self.offset_range[0] and
                                       abs(offsets[offset_exon_idx][1]) <= self.offset_range[1])
                if end_is_in_range:
                    entry_key = transcript_id + ".exon_" + \
                        str(offset_exon_idx) + ".end" + \
                        '.offset_' + str(offsets[offset_exon_idx][1])
                    for exon_number, exon in enumerate(self.reference_db.children(
                            reference_id, featuretype='exon', order_by='start')):
                        if exon_number == offset_exon_idx:
                            new_entry_key = (transcript_id, exon.end)
                            extracted_candidates[new_entry_key] = {
                                "transcript_id": transcript_id,
                                "strand": strand,
                                "exon_number": offset_exon_idx,
                                "location_type": "end",
                                "location": exon.end + offsets[offset_exon_idx][1],
                                "offset": offsets[offset_exon_idx][1],
                                "seq_id": exon.seqid
                            }
        output_manager.output_line({
            "line": f"Found {len(extracted_candidates)} candidates in the given offset range",
            "is_info": True
        })
        return extracted_candidates
