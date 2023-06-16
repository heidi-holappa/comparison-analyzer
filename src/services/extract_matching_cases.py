class MatchingCasesExtractor:

    def __init__(self, offset_results: dict, offset: int, reference_db):
        """
            Extract candidates matching the selected offset.

        Args:
            offset_results (dict): a dictionary containing the offset results for each transcript
            offset (int): the offset to match
            reference_db: a gffutils database with the reference annotation
        """
        self.offset_results = offset_results
        self.offset = offset
        self.reference_db = reference_db

    def extract_candidates_matching_selected_offset(self):
        extracted_candidates = {}
        for key, value in self.offset_results.items():
            if key[2] == '-':
                value.reverse()
            for i in range(1, len(value)-1):
                if abs(value[i][0]) == self.offset:
                    if value[i][0] > 0:
                        order_by = 'end'
                    else:
                        order_by = 'start'
                    for exon in self.reference_db.children(key[1], featuretype='exon', order_by=order_by):
                        if int(exon['exon_number'][0]) == i + 1:
                            extracted_candidates[(
                                key[0], key[1], key[2], i + 1, 'start')] = exon.start + value[i][0]
                            break
                elif abs(value[i][1]) == self.offset:
                    if value[i][1] > 0:
                        order_by = 'end'
                    else:
                        order_by = 'start'
                    for exon in self.reference_db.children(key[1], featuretype='exon', order_by=order_by):
                        if int(exon['exon_number'][0]) == i + 1:
                            extracted_candidates[(
                                key[0], key[1], key[2], i + 1, 'end')] = exon.end + value[i][1]
                            break
        return extracted_candidates
