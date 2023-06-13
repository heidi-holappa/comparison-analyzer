import pysam
from services.output_manager import default_output_manager as output_manager


class AlignmentParser:

    def __new__(cls):
        """
            Singleton pattern.
        """
        if not hasattr(cls, 'instance'):
            cls.instance = super(AlignmentParser, cls).__new__(cls)
        return cls.instance

    def __init__(self):
        """
            Parse a BAM-file and count the number of insertions and deletions at a given location.
            Follows the singleton pattern. Total case count for all processed BAM-files
            is stored in self.case_count.
        """
        self.case_count = {
            "insertion": 0,
            "deletion": 0,
        }

    def initialize_file(self, filename: str):
        # pylint: disable=no-member
        self.samfile = pysam.AlignmentFile(filename, "rb")

    def extract_location_from_cigar_string(self, cigar_tuples: list, reference_start: int, location: int):
        relative_position = location - reference_start
        alignment_position = 0
        ref_position = 0

        for idx, cigar_code in enumerate(cigar_tuples):
            if cigar_code[0] in [0, 2, 3, 7, 8]:
                ref_position += cigar_code[1]
            if ref_position < relative_position:
                alignment_position += cigar_code[1]
            else:
                return alignment_position + ref_position - relative_position
        return alignment_position + ref_position - relative_position

    def process_read(self, aligned_pairs: list, location: int):
        """
        Process a read and count the number of insertions and deletions at a given location.

        Args:
            location (int): given location
        """

        deletion_found = False
        insertion_found = False
        window_size = 8
        for element in aligned_pairs[location-window_size:location]:
            if not element[0] and not deletion_found:
                self.case_count["deletion"] += 1
                deletion_found = True
            if element[1] and not insertion_found:
                self.case_count["insertion"] += 1
                insertion_found = True
            if insertion_found and deletion_found:
                break

    def process_single_read(self, read, location: int):
        """
        Process a read and count the number of insertions and deletions at a given location.

        Args:
            location (int): given location
        """

        pairs = read.get_aligned_pairs()
        aligned = False
        deletion_found = False
        insertion_found = False
        for i in range(len(pairs) - 1):
            if pairs[i][0] and pairs[i][1] and pairs[i][1] < location:
                aligned = True

            if aligned and not pairs[i+1][1] == None and pairs[i+1][1] > location:
                # print(read.qname, pairs[i-5:i])
                for index in range(i-5, i):
                    if index < 0:
                        continue
                    if not pairs[index][0] and not deletion_found:
                        self.case_count["deletion"] += 1
                        deletion_found = True
                    if not pairs[index][1] and not insertion_found:
                        self.case_count["insertion"] += 1
                        insertion_found = True
                    if insertion_found and deletion_found:
                        break
                break

    def process_bam_file(self, reads_and_locations: dict):
        count = 0
        for read in self.samfile.fetch():
            if read.qname in reads_and_locations:
                count += 1
                if count % 1000 == 0:
                    output_manager.output_line({
                        "line": "Processed " + str(count) + " reads",
                        "end_line": "\r",
                        "is_info": True
                    })
                for location in reads_and_locations[read.qname]:
                    aligned_location = self.extract_location_from_cigar_string(
                        read.cigartuples,
                        read.reference_start,
                        location
                    )
                    print(aligned_location)
                    self.process_read(
                        read.get_aligned_pairs(), aligned_location)

    def execute(self, filename: str, reads_and_locations: dict):
        """
        Initialize AlignmentFile and execute the alignment parser.

        Args:
            filename (str): _description_
            location (int): _description_
        """
        output_manager.output_line({
            "line": "Processing file: " + filename,
            "is_info": True
        })
        self.initialize_file(filename)
        self.process_bam_file(reads_and_locations)


default_alignment_parser = AlignmentParser()
