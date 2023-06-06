from services.filter_bam_multifile import create_output_filename_dict, create_read_dict, filter_reads

class BamManager:

    def __init__(self, bam_path: str, tsv_path: str, transcript_set: set):
        self.bam_path = bam_path
        self.tsv_path = tsv_path
        self.transcript_set = transcript_set

    def execute(self):
        filename_dict = create_output_filename_dict(self.bam_path, self.transcript_set)
        read_dict = create_read_dict(self.bam_path, self.tsv_path)
        filter_reads(self.bam_path, filename_dict, read_dict)