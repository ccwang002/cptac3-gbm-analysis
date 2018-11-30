import gzip
from csv import DictReader

selected_samples = snakemake.params['selected_samples']

with gzip.open(snakemake.input[0], 'rt') as f:
    next(f)
    next(f)
    reader = DictReader(f, dialect='excel-tab')
    # make sure all the samples are there
    columns = ['Name', 'Description', *selected_samples]
    if not set(columns) <= set(reader.fieldnames):
        missing_samples = set(columns) - set(reader.fieldnames)
        raise ValueError(f"Cannot find {', '.join(missing_samples)} in the table.")

    with gzip.open(snakemake.output[0], 'wt') as outf:
        # Write header
        outf.write('\t'.join(columns))
        outf.write('\n')
        for row in reader:
            vals = [row[c] for c in columns]
            outf.write('\t'.join(vals))
            outf.write('\n')
