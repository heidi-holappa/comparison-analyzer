import pysam


class AlignmentParser:

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

    def process_read(self, location: int):
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

    def execute(self, filename: str, location: int):
        """
        Initialize AlignmentFile and execute the alignment parser.

        Args:
            filename (str): _description_
            location (int): _description_
        """
        self.initialize_file(filename)
        self.process_read(location)


default_alignment_parser = AlignmentParser()
