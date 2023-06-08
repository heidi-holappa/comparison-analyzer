import pysam


class CigarParser:

    def initialize_file(self, filename: str):
        # pylint: disable=no-member
        samfile = pysam.AlignmentFile(filename, "rb")

        return samfile

    def extract_cigar_symbol(self, samfile, coordinate):
        output = {}
        canonicals = {}
        for read in samfile.fetch(multiple_iterators=True):
            if read.reference_start > coordinate:
                continue
            modified_cigar = []
            ref_position = 0
            seq_position = 0
            cigar_tuples = read.cigartuples
            for i in range(len(cigar_tuples)):
                modified_cigar.append(
                    (ref_position, seq_position, cigar_tuples[i][0], cigar_tuples[i][1]))
                if cigar_tuples[i][0] in [0, 2, 3, 7, 8]:
                    ref_position += cigar_tuples[i][1]
                if cigar_tuples[i][0] in [0, 1, 4, 7, 8]:
                    seq_position += cigar_tuples[i][1]
            # print(read.query_name, read.query_alignment_start, modified_cigar)
            relative_position = coordinate - read.reference_start
            for i in range(len(modified_cigar)):
                if modified_cigar[i][0] <= relative_position and \
                        len(modified_cigar) > i+1 and modified_cigar[i+1][0] > relative_position:
                    local_point = relative_position - modified_cigar[i][0]
                    seq = read.query_sequence[modified_cigar[i][1] +
                                              local_point - 4:modified_cigar[i][1] + local_point]
                    if seq[2:] not in canonicals:
                        canonicals[seq[2:]] = 0
                    canonicals[seq[2:]] += 1

            output[(read.query_name, read.reference_start)] = modified_cigar

        print(len(output))
        print(canonicals)
