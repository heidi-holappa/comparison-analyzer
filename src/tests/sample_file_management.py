import os
from config import TEST_FILE_DIR

from services.db_initializer import init_databases


class SampleFileManagement:

    def __init__(self):
        self.gffcompare_gtf = os.path.join(
            TEST_FILE_DIR, "gffcompare.annotated.gtf")
        self.reference_gtf = os.path.join(
            TEST_FILE_DIR, "gencode.vM26.basic.annotation.extract.gtf")
        self.gffcompare_db_filename = os.path.join(
            TEST_FILE_DIR, "gffcompare.annotated-ca.db")
        self.reference_db_filename = os.path.join(
            TEST_FILE_DIR, "gencode.vM26.basic.annotation.extract-ca.db")
        self.bam_file = os.path.join(
            TEST_FILE_DIR, "Mouse.ONT.R9.4.sim.RE.no_gtf.transcript925.ch1.nnic.bam")
        self.gffcompare_db = None
        self.reference_db = None

    def initialize_test_files(self):
        self.gffcompare_db, self.reference_db = init_databases(
            self.gffcompare_gtf, self.reference_gtf, True)

    def databases_exist(self) -> bool:
        gffcompare_db = os.path.exists(self.gffcompare_db_filename)
        reference_db = os.path.exists(self.reference_db_filename)
        return gffcompare_db and reference_db

    def remove_test_files(self):
        if os.path.exists(self.gffcompare_db_filename):
            os.remove(self.gffcompare_db_filename)
        if os.path.exists(self.reference_db_filename):
            os.remove(self.reference_db_filename)


default_test_file_manager = SampleFileManagement()
