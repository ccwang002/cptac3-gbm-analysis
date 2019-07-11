"""Extract TinDaisy caller information as a matching column named callers"""
import argparse
import logging
from cyvcf2 import VCF

logger = logging.getLogger(__name__)


def read_caller_info(vcf_pth):
    vcf = VCF(vcf_pth)
    for i, mut in enumerate(vcf, 1):
        try:
            callers = mut.INFO['set']
            yield callers
        except ValueError as e:
            logger.exception(f'Caller information of {i}-th mutation is missing.')
            yield 'NA'
    logger.info(f'Read total {i:,d} mutations.')


def main(vcf_pth):
    print('callers')
    for callers in read_caller_info(vcf_pth):
        print(callers)


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
    parser.add_argument('vcf_pth', help="Path to the VCF")
    return parser


if __name__ == '__main__':
    parser = setup_cli()
    args = parser.parse_args()
    main(args.vcf_pth)
