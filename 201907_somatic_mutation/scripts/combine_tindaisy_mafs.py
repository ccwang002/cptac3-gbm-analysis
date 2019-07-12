import csv
import gzip
import itertools
import logging
from pathlib import Path
from maf_utils import MAF

logger = logging.getLogger(__name__)


def read_callers(callers_pth):
    with gzip.open(callers_pth, 'rt') as callers_f:
        assert next(callers_f).rstrip('\n') == 'callers'
        for line in callers_f:
            yield line.rstrip('\n')


def read_maf_content(maf_pth, callers_pth):
    maf_reader = MAF(maf_pth)
    callers_reader = read_callers(callers_pth)
    yield from zip(maf_reader, callers_reader)


def main(maf_pths, callers_pths, out_pth):
    logger.info(f'Combining {len(maf_pths):,d} MAFs.')
    # Make sure the same number of MAFs and callers files
    assert len(maf_pths) == len(callers_pths)
    # Read MAF header
    maf_reader = MAF(maf_pths[0])
    maf_cols = [*maf_reader.raw_columns, 'callers']

    # A chained iterator of all mutations
    mutation_iter = itertools.chain.from_iterable(
        map(read_maf_content, maf_pths, callers_pths)
    )
    with gzip.open(out_pth, 'wt', compresslevel=9) as outf:
        # Write header
        print('\t'.join(maf_cols), file=outf)
        for i, (mut, callers) in enumerate(mutation_iter, 1):
            if i % 10000 == 0:
                logger.info(f'... read {i:,d} mutation calls')
            print('\t'.join([*mut, callers]), file=outf)
        logger.info(f'... read total {i:,d} mutation calls')
    logger.info('Complete')


def setup_cli(snakemake):
    # Setup logging to file
    fh = logging.StreamHandler()
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(fh)
    log_fmt = '[%(asctime)s][%(levelname)-7s] %(message)s'
    log_formatter = logging.Formatter(log_fmt, '%Y-%m-%d %H:%M:%S')
    fh.setFormatter(log_formatter)


if __name__ == '__main__':
    setup_cli(snakemake)                    # noqa
    main(
        snakemake.input['maf'],             # noqa
        snakemake.input['callers'],         # noqa
        snakemake.output[0],                # noqa
    )
