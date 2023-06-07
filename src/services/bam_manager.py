import os
from services.filter_bam_multifile import create_output_filename_dict
from services.filter_bam_multifile import create_read_dict
from services.filter_bam_multifile import filter_reads


class BamManager:

    def __init__(self, bam_path: str, tsv_path: str, transcript_set: set):
        self.bam_path = bam_path
        self.tsv_path = tsv_path
        self.transcript_set = transcript_set

    def execute(self):
        output_filename_dict = create_output_filename_dict(
            self.bam_path,
            self.transcript_set,
            self.create_temporary_path()
        )
        read_dict = create_read_dict(
            output_filename_dict,
            self.tsv_path
        )
        filter_reads(
            self.bam_path,
            output_filename_dict,
            read_dict
        )

    def create_temporary_path(self):
        temporary_dir = "temporary_files"
        current_dir = os.getcwd()
        temporary_path = os.path.join(current_dir, temporary_dir)
        if not os.path.exists(temporary_path):
            os.makedirs(temporary_path)
        return temporary_path
