import gffutils
import argparse

parser = argparse.ArgumentParser(description='comAna: a tool for comparing annotations')
parser.add_argument('-i', '--input', help='gtf-file to be imported into the database', required=True)


arguments = parser.parse_args()
