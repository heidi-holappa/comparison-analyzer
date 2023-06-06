from pyfaidx import Fasta

class FastaExtractor:

    def __init__(self, fasta_path: str):
        self.fasta = Fasta(fasta_path)

    def extract_characters_at_given_coordinates(self, coordinates: tuple):
        chromosome, start, end = coordinates
        return self.fasta[chromosome][start:end]

    


# genes = Fasta('/home/holaphei/koulutyot/lv-2022-2023/BI-summer-trainee/data/sprint-1-get-familiar-with-isoquant/reference/mouse/GRCm39.primary_assembly.genome.fa')
# print(genes.keys())
# s = genes['chr6'][87866109-4:87866109+4]
# print(s)
# ENSMUST00000068755.14 