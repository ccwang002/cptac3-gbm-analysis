"""Extract and transpose methylation data.

It converts the given CSV file:

              probe_1   probe_2   ...
    sample_A
    sample_B

to be transposed:
              sample_A  sample_B
    probe_1
    probe_2
    ...
"""
import argparse
import csv
import gzip
from pathlib import Path
import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def convert_to_float(x):
    """Convert an element to float that is aware of NA."""
    if x == 'NA':
        return float('nan')
    else:
        return float(x)


def transpose_meth_table(reader):
    probe_ids = next(reader)[1:]
    df_dict = {'probe_id': probe_ids}
    for i, row in enumerate(reader, 1):
        sample_id = row[0]
        logger.info(f'Converting {i}-th sample: {sample_id}')
        # Alternatively, use:
        # arr = pd.to_numeric(one_sample[1:], errors='coerce')
        beta_arr = np.fromiter(map(convert_to_float, row[1:]), dtype=np.float64)
        df_dict[sample_id] = beta_arr

    # Convert it to the dataframe
    df = pd.DataFrame.from_dict(df_dict)
    return df


def main(csv_pth, out_pth):
    with gzip.open(csv_pth, 'rt') as f:
        reader = csv.reader(f, dialect='excel')
        df = transpose_meth_table(reader)
        df.to_feather(out_pth)


if __name__ == '__main__':
    # Setup console logging
    console = logging.StreamHandler()
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(console)
    log_formatter = logging.Formatter('%(message)s')
    console.setFormatter(log_formatter)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__
    )
    parser.add_argument('csv', help="Path to Methylation CSV file")
    parser.add_argument('out', default='.', help="Output folder")
    args = parser.parse_args()
    main(args.csv, args.out)
