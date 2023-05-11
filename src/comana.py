import os
import gffutils
import argparse

parser = argparse.ArgumentParser(description='comAna: a tool for comparing annotations')
parser.add_argument('-i', '--input', help='gtf-file to be imported into the database', required=True)


arguments = parser.parse_args()


db_path = os.path.dirname(arguments.input)
db_name = os.path.basename(arguments.input)[:-4] + '-ca.db'
print(db_name)
print(db_path)


# Sample file path
# /home/holaphei/koulutyot/lv-2022-2023/BI-summer-trainee/isoquant-summer-trainee/ca-task-2/comparison-analyzer/sample-data/nanopore-compare-2.annotated.gtf