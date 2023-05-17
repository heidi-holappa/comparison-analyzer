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
| c | class-code | specify one or several class codes for which to generate offset data. | 

**Example**

The following instruction computes analysis for the files given with arguments `i` and `r`. If sqlite-databases exist for the given files, those are used. Argument `s` specifies that overview data is to be genegated and argument `-c` specifies that offsets should be calculated for classcodes `i` and `x`.

```
python3 src/compana.py -i <gffcompare-file> -r <reference-file> -s -c i x
```

The following instruction takes files given with arguments `i` and `r`. Argument `f` specifies that new databases are to be created using `gffutils` from the given files. Any existing databases will overwritten.

```
python3 src/compana.py -i <gffcompare-file> -r <reference-file> -f
```