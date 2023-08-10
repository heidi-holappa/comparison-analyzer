import os
import sys
from datetime import datetime
from dotenv import load_dotenv

dirname = os.path.dirname(__file__)

try:
    load_dotenv(dotenv_path=os.path.join(dirname, "..", ".env"))
except FileNotFoundError:
    print("Warning! .env file not found.")

temporary_dir_path = os.getenv("TEMPORARY_DIR") or "temporary_files"
offset_log = os.getenv("OFFSET_LOG_FILE") or "debug_offsets.log"
test_file_directory = os.getenv(
    "TEST_FILE_DIRECTORY") or "src/tests/files"

saved_results = os.getenv("SAVED_RESULTS_DIR") or "saved_results"
log_dir = os.getenv("LOG_FILE_DIR") or "logs"
fasta_overview = os.getenv("FASTA_OVERVIEW_FILE") or "fasta_overview.md"
cigar_results = os.getenv("CIGAR_RESULTS_LOG") or "cigar_results.md"


LOG_DIR = os.path.join(dirname, "..", log_dir)

SAVE_DIR = os.path.join(dirname, "..", "saved_results")

if len(sys.argv) > 1 and sys.argv[1] == "-j":
    log_extension = os.path.splitext(sys.argv[2])[0]
    if '-n' in sys.argv:
        log_extension += "-no_canonicals"
else:
    log_extension = ""

LOG_FILE_DIR = os.path.join(dirname, "..", log_dir,
                            log_extension + "-" + datetime.now().strftime('%Y-%m-%d_%H-%M'))
TEMPORARY_DIR = os.path.join(dirname, "..", temporary_dir_path)
OFFSET_LOG = os.path.join(LOG_FILE_DIR, offset_log)
TEST_FILE_DIR = os.path.join(dirname, "..", test_file_directory)
FASTA_OVERVIEW_FILE = os.path.join(LOG_FILE_DIR, fasta_overview)
CIGAR_RESULTS_LOG = os.path.join(LOG_FILE_DIR, cigar_results)
TITLE_FILE_LENGTH = os.getenv("TITLE_FILE_LENGTH") or 50
DEFAULT_WINDOW_SIZE = os.getenv("DEFAULT_WINDOW_SIZE") or 8
CREATE_IMG_N_TRESHOLD = os.getenv("CREATE_IMG_N_TRESHOLD") or 100

if not os.path.exists(LOG_DIR):
    os.mkdir(LOG_DIR)
if not os.path.exists(LOG_FILE_DIR):
    os.mkdir(LOG_FILE_DIR)
if not os.path.exists(TEMPORARY_DIR):
    os.mkdir(TEMPORARY_DIR)
if not os.path.exists(SAVE_DIR):
    os.mkdir(SAVE_DIR)
