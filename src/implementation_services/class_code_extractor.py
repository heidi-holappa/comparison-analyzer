
def extract_class_codes_with_transcripts(db_gffcompare):
    class_code_dict = {}
    for transcript in db_gffcompare.features_of_type('transcript'):
        if 'class_code' in transcript.attributes:
            if transcript.id not in class_code_dict:
                class_code_dict[transcript.id] = transcript.attributes['class_code'][0]
    return class_code_dict
