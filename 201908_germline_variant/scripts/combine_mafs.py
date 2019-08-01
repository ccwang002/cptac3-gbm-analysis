import argparse
import csv
import gzip
import itertools
import logging
from pathlib import Path
from maf_utils import MAF

logger = logging.getLogger(__name__)


def read_maf_content(maf_pth):
    logger.info(f'... opening new MAF at {maf_pth}')
    maf_reader = MAF(maf_pth)
    yield from maf_reader


def main(maf_folder, out_pth):
    maf_pths = sorted(Path(maf_folder).glob('*.maf.gz'))
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


def setup_cli():
    # Setup console logging
    console = logging.StreamHandler()
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(console)
    log_fmt = '[%(asctime)s][%(levelname)-7s] %(message)s'
    log_formatter = logging.Formatter(log_fmt, '%Y-%m-%d %H:%M:%S')
    console.setFormatter(log_formatter)

    parser = argparse.ArgumentParser()
    parser.add_argument('maf_folder', help="Path to folder contains all MAFs")
    parser.add_argument('out_maf_gz', help="Path to output maf.gz")
    return parser


if __name__ == '__main__':
    parser = setup_cli()
    args = parser.parse_args()
    main(args.maf_folder, args.out_maf_gz)
