from config import DEFAULT_WINDOW_SIZE
from services.output_manager import default_output_manager as output_manager


class CaseExtractor:

    def __init__(self):
        """
            Extract all intron site locations from isoquant gtf-db.

        Args:
            isoquant_db: reference to a gffutils database with the isoquant annotation
        """
        pass

    def extract_intron_site_locations(self, isoquant_db, window_size: int = int(DEFAULT_WINDOW_SIZE)) -> dict:
        output_manager.output_line({
            "line": "IMPLEMENTATION: EXTRACTING CASES",
            "is_title": True
        })
        output_manager.output_line({
            "line": f"Extracting intron site locations from the IsoQuant gtf-db",
            "is_info": True
        })
        extracted_cases = {}

        for transcript in isoquant_db.features_of_type('transcript'):
            for exon in isoquant_db.children(transcript, featuretype='exon', order_by='start'):
                transcript_id = str(transcript.id)
                case_information = {
                    'exon_start': {
                        'case_key': (transcript_id, str(exon.start)),
                        'location': exon.start,
                        'location_type': 'start'
                    },
                    'exon_end': {
                        'case_key': (transcript_id, str(exon.end)),
                        'location': exon.end,
                        'location_type': 'end'
                    }
                }
                for case in case_information.values():
                    extracted_cases[case['case_key']] = {
                        'transcript_id': transcript_id,
                        'strand': transcript.strand,
                        'location_type': case['location_type'],
                        "location": case['location'],
                        "seq_id": exon.seqid,
                        "extracted_information": {
                            "left": {
                                "insertions": {},
                                "deletions": {},
                                "ins_pos_distr": [0] * window_size,
                                "del_pos_distr": [0] * window_size,
                                "closest_canonical": "",
                                "error_detected": False

                            },
                            "right": {
                                "insertions": {},
                                "deletions": {},
                                "ins_pos_distr": [0] * window_size,
                                "del_pos_distr": [0] * window_size,
                                "closest_canonical": "",
                                "error_detected": False
                            }
                        }
                    }

        output_manager.output_line({
            "line": f"Total of {len(extracted_cases)} intron sites extracted from isoquant gtf-db",
            "is_info": True
        })
        return extracted_cases

    def extract_transcripts(self, isoquant_db) -> set:
        transcripts = set()

        for transcript in isoquant_db.features_of_type('transcript'):
            transcripts.add(str(transcript.id))

        return transcripts
