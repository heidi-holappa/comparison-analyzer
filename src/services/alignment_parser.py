import pysam


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

    def process_single_read(self, location: int):
        """
        Process a read and count the number of insertions and deletions at a given location.

        Args:
            location (int): given location
        """
        for read in self.samfile.fetch():
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
                    break

    def process_bam_file(self, reads_and_locations: dict):
        for read in self.samfile.fetch():
            if read.qname in reads_and_locations:
                for location in reads_and_locations[read.qname]:
                    self.process_single_read(location)

    def execute(self, filename: str, reads_and_locations: dict):
        """
        Initialize AlignmentFile and execute the alignment parser.

        Args:
            filename (str): _description_
            location (int): _description_
        """
        self.initialize_file(filename)
        self.process_bam_file(reads_and_locations)


default_alignment_parser = AlignmentParser()