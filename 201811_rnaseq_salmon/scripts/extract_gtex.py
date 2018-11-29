import gzip
from csv import DictReader

with gzip.open(snakemake.input[0], 'rt') as f:
    next(f)
    next(f)
    reader = DictReader(f, dialect='excel-tab')
    # make sure all the samples are there

    for row in reader:

