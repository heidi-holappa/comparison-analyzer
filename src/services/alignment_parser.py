import pysam
from services.output_manager import default_output_manager as output_manager
from config import DEFAULT_WINDOW_SIZE


class AlignmentParser:

    def __new__(cls):
        """
            Singleton pattern.
        """
        if not hasattr(cls, 'instance'):
            cls.instance = super(AlignmentParser, cls).__new__(cls)
        return cls.instance

    # TODO: remove reads and transcripts. It is for debugging purposes.
    def __init__(self):
        """
            Parse a BAM-file and count the number of insertions and deletions at a given location.
            Follows the singleton pattern. Total case count for all processed BAM-files
            is stored in self.case_count.
        """

        self.reads_and_transcripts = {}
        self.case_count = {
            "insertions": {},
            "deletions": {},
        }
        # TODO: add to configuration or to user arguments (window size)
        self.window_size = int(DEFAULT_WINDOW_SIZE)

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
            if ref_position <= relative_position:
                alignment_position += cigar_code[1]
            else:
                return alignment_position + (cigar_code[1] - (ref_position - relative_position))

        return -1

    def process_read(self, aligned_pairs: list, location: int, loc_type: str):
        """
        Process a read and count the number of insertions and deletions at a given location.

        Args:
            aligned_pairs (list): list of tuples of aligned pairs
            location (int): given location
            type (str): type of location (start or end)
        """

        deletions = 0
        insertions = 0
        debug_list = []

        if loc_type == "end":
            debug_list = aligned_pairs[location - self.window_size:location]
            for element in aligned_pairs[location - self.window_size:location]:
                if not element[0]:
                    deletions += 1
                if not element[1]:
                    insertions += 1
        elif loc_type == "start":
            debug_list = aligned_pairs[location:location + self.window_size]
            for element in aligned_pairs[location:location + self.window_size]:
                if not element[0]:
                    deletions += 1
                if not element[1]:
                    insertions += 1

        if deletions:
            if deletions not in self.case_count["deletions"]:
                self.case_count["deletions"][deletions] = 0
            self.case_count["deletions"][deletions] += 1
        if insertions:
            if insertions not in self.case_count["insertions"]:
                self.case_count["insertions"][insertions] = 0
            self.case_count["insertions"][insertions] += 1

        if deletions > 8 or insertions >= 8:
            return True, debug_list
        return False, debug_list

    def process_bam_file(self, reads_and_locations: dict):
        count = 0
        errors = []
        for read in self.samfile.fetch():
            if read.is_supplementary:
                continue
            if read.query_name in reads_and_locations:
                count += 1
                if count % 1000 == 0:
                    output_manager.output_line({
                        "line": "Processed " + str(count) + " reads",
                        "end_line": "\r",
                        "is_info": True
                    })
                for location, type in reads_and_locations[read.query_name]:
                    location = location - 1

                    if read.reference_start > location or read.reference_end < location:
                        continue

                    if not read.cigartuples:
                        continue

                    aligned_location = self.extract_location_from_cigar_string(
                        read.cigartuples,
                        read.reference_start,
                        location
                    )

                    response, debug_list = self.process_read(
                        read.get_aligned_pairs(),
                        aligned_location,
                        type)
                    if response:
                        errors.append(
                            f"{read.query_name}\t{self.reads_and_transcripts[read.query_name]}\t{location}\t{aligned_location}\t{type}\t{read.reference_start}\t{read.reference_end}\t{debug_list}\n")
        with open("alignment_errors.txt", "w") as file:
            file.write(
                "qname\ttranscripts\tlocation\talign_location\ttype\tread.reference_start\tread.reference_end\tlist of alignments\n")
            file.writelines(errors)

        output_manager.output_line({
            "line": "\nFinished",
            "is_info": True
        })

    # TODO: remove transcripts and reads. It is here for debugging
    def execute(self, filename: str, reads_and_locations: dict, transcripts_and_reads: dict):
        """
        Initialize AlignmentFile and execute the alignment parser.

        Args:
            filename (str): _description_
            location (int): _description_
        """

        for key, value in transcripts_and_reads.items():
            for read in value:
                if read not in self.reads_and_transcripts:
                    self.reads_and_transcripts[read] = set()
                self.reads_and_transcripts[read].add(key)

        self.initialize_file(filename)
        self.process_bam_file(reads_and_locations)


default_alignment_parser = AlignmentParser()
