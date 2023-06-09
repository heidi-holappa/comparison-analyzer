import os
from dotenv import load_dotenv

dirname = os.path.dirname(__file__)

try:
    load_dotenv(dotenv_path=os.path.join(dirname, "..", ".env"))
except FileNotFoundError:
    print("Warning! .env file not found.")

temporary_dir_path = os.getenv("TEMPORARY_DIR") or "temporary_files"
offset_log = os.getenv("OFFSET_LOG") or "offset_log.txt"
test_file_directory = os.getenv(
    "TEST_FILE_DIRECTORY") or "src/tests/files"

log_dir = os.getenv("LOG_FILE_DIR") or "logs"
fasta_overview = os.getenv("FASTA_OVERVIEW_FILE") or "fasta_overview.md"
cigar_results = os.getenv("CIGAR_RESULTS_LOG") or "cigar_results.md"
LOG_FILE_DIR = os.path.join(dirname, "..", log_dir)
if not os.path.exists(LOG_FILE_DIR):
    os.mkdir(LOG_FILE_DIR)

TEMPORARY_DIR = os.path.join(dirname, "..", temporary_dir_path)
OFFSET_LOG = os.path.join(LOG_FILE_DIR, offset_log)
TEST_FILE_DIR = os.path.join(dirname, "..", test_file_directory)
FASTA_OVERVIEW_FILE = os.path.join(LOG_FILE_DIR, fasta_overview)
CIGAR_RESULTS_LOG = os.path.join(LOG_FILE_DIR, cigar_results)
TITLE_FILE_LENGTH = os.getenv("TITLE_FILE_LENGTH") or 50
