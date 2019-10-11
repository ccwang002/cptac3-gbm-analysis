"""Extract TCGA MC3 mutation calls LGG only."""
import argparse
import csv
import gzip
import logging
from maf_utils import MAF

logger = logging.getLogger(__name__)


class MC3MAF(MAF):
    def make_columns(self, raw_columns):
        # Strand is duplicated
        renamed_cols = []
        for col in raw_columns:
            # Follows GDC MAF standard to rename VEP's STRAND column to be TRANSCRIPT_STRAND
            if col == 'STRAND':
                renamed_cols.append('TRANSCRIPT_STRAND')
            else:
                renamed_cols.append(col)
        return renamed_cols


def get_all_tcga_lgg_cases(tsv_pth):
    with gzip.open(tsv_pth, 'rt') as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        all_cases = set(row['submitter_id'] for row in reader)
    return all_cases


def get_all_tcga_lgg_tss_codes(tsv_pth):
    # TSS = Tissue Source Site
    with gzip.open(tsv_pth, 'rt') as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        all_lgg_tss_codes = set(
            row['TSS Code'] for row in reader
            if row['Study Name'] == 'Brain Lower Grade Glioma'
        )
    return all_lgg_tss_codes


def main(maf_pth, case_pth, tss_pth):
    logger.info(
        f'Read MC3 MAF from {maf_pth}, all lgg cases from {case_pth}, TCGA TSS code from {tss_pth}'
    )
    maf = MC3MAF(maf_pth)
    all_tcga_lgg_cases = get_all_tcga_lgg_cases(case_pth)
    all_tcga_lgg_tss_codes = get_all_tcga_lgg_tss_codes(tss_pth)

    seen_lgg_cases = set()

    # Print out the header
    print('\t'.join(maf.columns))
    for i, mut in enumerate(maf, 1):
        if i % 50000 == 0:
            logger.info(f'... processed {i:,d} mutations.')

        case_id = mut.Tumor_Sample_Barcode[:len('TCGA-XX-XXXX')]
        tss = case_id[len('TCGA-'):len('TCGA-XX')]

        is_lgg_tss = tss in all_tcga_lgg_tss_codes
        is_lgg_case = case_id in all_tcga_lgg_cases
        if is_lgg_tss:
            if not is_lgg_case:
                logger.error(f'{i}-th mutation has lgg TSS {tss_id} but not a LGG case {case_id}')
            seen_lgg_cases.add(case_id)
            # Export the mutation
            print('\t'.join(mut))
        elif is_lgg_case:
            logger.error(f'{i}-th mutation has no lgg TSS {tss_id} but is a LGG case {case_id}')

    logger.info(f'All mutations are processed. '
                f'Detected {len(seen_lgg_cases)} LGG cases with MC3 mutations.')


def setup_cli():
    # Setup console logging
    console = logging.StreamHandler()
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(console)
    log_fmt = '[%(asctime)s][%(levelname)-7s] %(message)s'
    log_formatter = logging.Formatter(log_fmt, '%Y-%m-%d %H:%M:%S')
    console.setFormatter(log_formatter)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('maf', help="Path to the MC3 MAF")
    parser.add_argument('case_tsv', help="Path to all TCGA lgg cases")
    parser.add_argument('tss_tsv', help="Path to TCGA TSS code")
    return parser


if __name__ == '__main__':
    parser = setup_cli()
    args = parser.parse_args()
    main(args.maf, args.case_tsv, args.tss_tsv)
