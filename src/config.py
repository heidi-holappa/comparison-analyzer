import os
from dotenv import load_dotenv

dirname = os.path.dirname(__file__)

try:
    load_dotenv(dotenv_path=os.path.join(dirname, "..", ".env"))
except FileNotFoundError:
    print("Warning! .env file not found.")

temporary_dir_path = os.getenv("TEMPORARY_DIR") or "temporary_files"
offset_log_path = os.getenv("OFFSET_LOG") or "offset_log.txt"

TEMPORARY_DIR = os.path.join(dirname, "..", temporary_dir_path)
OFFSET_LOG = os.path.join(dirname, "..", offset_log_path)
