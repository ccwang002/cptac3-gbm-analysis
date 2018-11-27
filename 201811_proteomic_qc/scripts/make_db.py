import argparse
from collections import namedtuple
from typing import NamedTuple
import logging
from pathlib import Path
import re
from sqlalchemy import (
    create_engine, event,
    MetaData, Table, Column, Integer, Text, Index,
    UniqueConstraint
)
from sqlalchemy.engine import Engine
from cyvcf2 import VCF

logger = logging.getLogger(__name__)
BATCH_SIZE = 5000

GROUP_GDC_RAW_CALLS_SQL = """\
CREATE TABLE gdc_raw AS
SELECT
    sample, chromosome, start, end, ref_allele, alt_allele,
    avg(t_depth) AS t_depth,
    avg(t_ref_count) AS t_ref_count,
    avg(t_alt_count) AS t_alt_count,
    avg(n_depth) AS n_depth,
    avg(n_ref_count) AS n_ref_count,
    avg(n_alt_count) AS n_alt_count,
    avg(t_alt_count) / avg(t_depth) AS t_allele_freq,
    avg(n_alt_count) / avg(n_depth) AS n_allele_freq,
    group_concat(caller, '|') AS callers
FROM gdc_raw_per_caller
WHERE filter IS NULL
GROUP BY sample, chromosome, start, end, ref_allele, alt_allele
"""


class RawMutationCall(NamedTuple):
    sample: str
    chromosome: str
    start: int
    end: int
    ref_allele: str
    alt_allele: str
    filter: str
    caller: str
    t_depth: int
    t_ref_count: int
    t_alt_count: int
    n_depth: int
    n_ref_count: int
    n_alt_count: int


class AnnotatedMutationCall(NamedTuple):
    sample: str
    chromosome: str
    start: int
    end: int
    ref_allele: str
    alt_allele: str
    filter: str
    caller: str
    t_depth: int
    t_ref_count: int
    t_alt_count: int
    n_depth: int
    n_ref_count: int
    n_alt_count: int
    symbol: str
    transcript: str
    tls: int
    biotype: str
    consequence: str
    hgvs_c: str
    hgvs_p: str


def cal_read_counts(mut, caller):
    """Extract the read count from the VCF."""
    if caller == 'muse' or caller == 'mutect2':
        n_depth, t_depth = mut.gt_depths
        n_ref_count, t_ref_count = mut.gt_ref_depths
        n_alt_count, t_alt_count = mut.gt_alt_depths
    elif caller == 'varscan2':
        n_depth, t_depth = mut.format('DP')[:, 0]
        n_ref_count, t_ref_count = mut.format('RD')[:, 0]
        n_alt_count, t_alt_count = mut.format('AD')[:, 0]
    elif caller == 'somaticsniper':
        n_depth, t_depth = mut.format('DP')[:, 0]
        n_ref_count, t_ref_count = mut.format('DP4')[:, :2].sum(axis=1)
        n_alt_count, t_alt_count = mut.format('DP4')[:, 2:].sum(axis=1)
    return t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count


def read_raw_vcf(vcf_p):
    """Read raw VCF"""
    sample, caller, *_ = Path(vcf_p).name.split('.')
    for mut in VCF(str(vcf_p)):
        # Convert chromosome
        if mut.CHROM == 'chrM':
            chrom = 'MT'
        else:
            chrom = mut.CHROM[3:]

        read_counts = list(map(int, cal_read_counts(mut, caller)))

        yield RawMutationCall(
            sample,
            chrom, mut.start, mut.end,
            mut.REF, mut.ALT[0],
            mut.FILTER,
            caller,
            *read_counts
        )



def create_cols(namedtuple_cls):
    """Create SQLAlchemy columns based on a named tuple class."""
    columns = []
    for field, f_type in namedtuple_cls._field_types.items():
        if f_type is int:
            col = Column(field, Integer())
        else:
            col = Column(field, Text())
        columns.append(col)
    return columns


def define_db_schema(metadata):
    Table(
        'gdc_raw_per_caller', metadata,
        *create_cols(RawMutationCall)
    )


def read_vcfs(vcf_pths):
    for vcf_p in vcf_pths:
        logger.info(f'... reading {vcf_p}')
        yield from read_raw_vcf(vcf_p)


def load_vcfs_to_db(conn, metadata, vcf_pths):
    ins = metadata.tables['gdc_raw_per_caller'].insert()
    # Insert mutation calls by batches
    ins_batch = []
    for i, mut in enumerate(read_vcfs(vcf_pths), 1):
        # Add new record into batch
        ins_batch.append(mut._asdict())

        # Commit the batch when a batch is full
        if len(ins_batch) >= BATCH_SIZE:
            with conn.begin():
                conn.execute(ins, ins_batch)
            ins_batch = []

        if i % 50000 == 0:
            logger.info(f'... processed {i:,d} mutation calls')

    # Load the last batch
    if ins_batch:
        with conn.begin():
            conn.execute(ins, ins_batch)
    logger.info(f'Loaded total {i:,d} mutation calls')


def main(db_pth, gdc_raw_vcfs):
    logger.info(f'Load mutation calls to SQLite at {db_pth}')
    logger.info(f'... loading {len(gdc_raw_vcfs)} VCFs')
    metadata = MetaData()
    db_engine = create_engine('sqlite:///' + db_pth)

    define_db_schema(metadata)
    metadata.create_all(db_engine)

    conn = db_engine.connect()
    load_vcfs_to_db(conn, metadata, gdc_raw_vcfs)

    logger.info(f'Group all GDC raw calls')
    conn.execute(GROUP_GDC_RAW_CALLS_SQL)

    logger.info('Complete')


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA cache_size=-8000000")
    cursor.execute("PRAGMA temp_store=MEMORY")
    cursor.execute("PRAGMA journal_mode=MEMORY")
    cursor.close()


def setup_cli():
    # Setup console logging
    console = logging.StreamHandler()
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(console)
    log_fmt = '[%(asctime)s][%(levelname)-7s] %(message)s'
    log_formatter = logging.Formatter(log_fmt, '%Y-%m-%d %H:%M:%S')
    console.setFormatter(log_formatter)


def pick_effect(mut, all_effects):
    """Pick the most representative effect of a mutation."""
    pass


def read_annotated_vcf(vcf_p):
    """Read annotated VCF"""
    sample, caller, *_ = Path(vcf_p).name.split('.')
    vcf = VCF(str(vcf_p))

    # Define a class to store the annotation record
    csq_header = next(header for header in vcf.header_iter() if header['ID'] == 'CSQ')
    csq_format = re.search(r'Format: ([\w\|]+)[\'"]', csq_header['Description']).group(1)
    CSQ = namedtupled('CSQ', csq_format.split('|'))

    for mut in vcf:
        # Convert chromosome
        if mut.CHROM == 'chrM':
            chrom = 'MT'
        else:
            chrom = mut.CHROM[3:]

        # Obtain read counts
        read_counts = list(map(int, cal_read_counts(mut, caller)))

        # Convert all the reported effects into CSQ objects
        all_effects = [
            CSQ(*eff_str.replace('&', ',').split('|'))
            for eff_str in mut.INFO['CSQ'].split(',')
        ]

        # Pick one effect to describe the mutation
        eff = pick_effect(mut, all_effects)

        yield RawMutationCall(
            sample,
            chrom, mut.start, mut.end,
            mut.REF, mut.ALT[0],
            mut.FILTER,
            caller,
            *read_counts,
            *eff
        )


if __name__ == '__main__':
    setup_cli()
    main(snakemake.output['db'], snakemake.input['raw_vcfs'])
