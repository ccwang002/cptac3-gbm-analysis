import csv
import gzip
import itertools
import logging
from pathlib import Path
from maf_utils import MAF

logger = logging.getLogger(__name__)


def read_maf_content(maf_pth):
    maf_reader = MAF(maf_pth)
    yield from maf_reader


def main(maf_pths, out_pth):
    logger.info(f'Combining {len(maf_pths):,d} MAFs.')
    # Read MAF header
    maf_reader = MAF(maf_pths[0])
    maf_cols = maf_reader.raw_columns

    # A chained iterator of all mutations
    mutation_iter = itertools.chain.from_iterable(
        map(read_maf_content, maf_pths)
    )
    with gzip.open(out_pth, 'wt', compresslevel=9) as outf:
        # Write header
        print('\t'.join(maf_cols), file=outf)
        for i, mut in enumerate(mutation_iter, 1):
            if i % 50000 == 0:
                logger.info(f'... read {i:,d} variant calls')
            print('\t'.join(mut), file=outf)
        logger.info(f'... read total {i:,d} variant calls')
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
        snakemake.input['all_mafs'],        # noqa
        snakemake.output[0],                # noqa
    )

