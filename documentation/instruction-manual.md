# Instruction manual

## Installation
Clone the project from GitHub. Create and launch a virtual environment for python 

```
python3 -m venv venv
source venv/bin/activate
```

and install dependencies
```
pip install -r requirements.txt
```

## Running application

The application provides the following arguments

| Short | Long | Description |
| --- | --- | --- |
| i | input | provide full path for `gtf` file produced by `gffcompare` |
| r | reference | provide full path for reference `gtf` file used to create the `gffcompare` `gtf`-file |
| f | force | force re-creation of sqlite3 database. By default a new database is not created if one already exists to improve efficiency. |
| s | stats | output statistics on class codes | 