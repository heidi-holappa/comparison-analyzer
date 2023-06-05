from pyfaidx import Fasta
genes = Fasta('/home/holaphei/koulutyot/lv-2022-2023/BI-summer-trainee/data/sprint-1-get-familiar-with-isoquant/reference/mouse/GRCm39.primary_assembly.genome.fa')
print(genes.keys())
s = genes['chr6'][87866109-4:87866109+4]
print(s)
# ENSMUST00000068755.14 