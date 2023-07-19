from services.output_manager import default_output_manager as output_manager


class CaseExtractor:

    def __init__(self):
        """
            Extract all intron site locations from isoquant gtf-db.

        Args:
            isoquant_db: reference to a gffutils database with the isoquant annotation
        """
        pass

    def extract_intron_site_locations(self, isoquant_db) -> dict:
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
                exon_start_entry_key = transcript_id + str(exon.start)
                exon_end_entry_key = transcript_id + str(exon.end)
                extracted_cases[exon_start_entry_key] = {
                    'transcript_id': transcript_id,
                    'strand': transcript.strand,
                    'location_type': "start",
                    "location": exon.start
                }
                extracted_cases[exon_end_entry_key] = {
                    'transcript_id': transcript_id,
                    'strand': transcript.strand,
                    'location_type': "end",
                    "location": exon.end
                }

        output_manager.output_line({
            "line": f"Total of {len(extracted_cases)} intron sites extracted from isoquant gtf-db",
            "is_info": True
        })
        return extracted_cases

    def extract_transcripts(self, isoquant_db) -> set:
        transcripts = set()

        for transcript in isoquant_db.features_of_type('transcript'):
            transcripts.add(transcript)

        return transcripts
