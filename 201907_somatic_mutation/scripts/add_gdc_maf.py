import logging
from pathlib import Path
from typing import NamedTuple
from sqlalchemy import (
    create_engine, event,
    MetaData, Table, Column, Integer, Text,
)
from sqlalchemy.engine import Engine
from cyvcf2 import VCF
from maf_utils import GDCMAF

logger = logging.getLogger(__name__)
BATCH_SIZE = 5000

GROUP_GDC_ANNOTATED_CALLS_SQL = """\
CREATE TABLE gdc AS
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
    group_concat(caller, '|') AS callers,
    group_concat("filter", '|') AS "filter",
    variant_type, variant_classification,
    symbol, gene, transcript_id, tsl, biotype, consequence, canonical,
    hgvsc, hgvsp, hgvsp_short,
    all_effects
FROM gdc_per_caller
WHERE filter IN ('PASS', 'panel_of_normals')
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
    return [int(x) for x in [t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count]]


def read_all_raw_vcfs(vcf_pths):
    """Read raw GDC VCF"""
    for vcf_p in vcf_pths:
        logger.info(f'... reading {vcf_p}')
        sample, caller, *_ = Path(vcf_p).name.split('.')
        for mut in VCF(str(vcf_p)):
            # Convert chromosome
            if mut.CHROM == 'chrM':
                chrom = 'MT'
            else:
                chrom = mut.CHROM[3:]
            # Get read counts
            read_counts = cal_read_counts(mut, caller)
            yield RawMutationCall(
                sample, chrom, mut.start, mut.end,
                mut.REF, mut.ALT[0], mut.FILTER, caller,
                *read_counts
            )


def load_vcfs_to_db(conn, metadata, vcf_pths):
    ins = metadata.tables['gdc_raw_per_caller'].insert()
    # Insert mutation calls by batches
    ins_batch = []
    for i, mut in enumerate(read_all_raw_vcfs(vcf_pths), 1):
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
    variant_type: str
    variant_classification: str
    symbol: str
    gene: str
    transcript_id: str
    tsl: str
    biotype: str
    consequence: str
    canonical: str
    hgvsc: str
    hgvsp: str
    hgvsp_short: str
    all_effects: str


def read_all_mafs(maf_pths):
    def to_int(x):
        if x == '':
            return 0
        else:
            return int(x)

    for maf_p in maf_pths:
        logger.info(f'... reading {maf_p}')
        sample, caller, *_ = Path(maf_p).name.split('.')
        maf = GDCMAF(maf_p, caller=caller)
        for m in maf:
            yield AnnotatedMutationCall(
                sample=sample, chromosome=m.chromosome,
                start=int(m.start_position), end=int(m.end_position),
                ref_allele=m.reference_allele, alt_allele=m.tumor_seq_allele2,
                filter=m.filter, caller=m.caller,
                t_depth=to_int(m.t_depth), t_ref_count=to_int(m.t_ref_count),
                t_alt_count=to_int(m.t_alt_count),
                n_depth=to_int(m.n_depth), n_ref_count=to_int(m.n_ref_count),
                n_alt_count=to_int(m.n_alt_count),
                variant_type=m.variant_type, variant_classification=m.variant_classification,
                symbol=m.symbol, gene=m.gene, transcript_id=m.transcript_id, tsl=m.tsl,
                biotype=m.biotype, consequence=m.consequence, canonical=m.canonical,
                hgvsc=m.hgvsc, hgvsp=m.hgvsp, hgvsp_short=m.hgvsp_short,
                all_effects=m.all_effects
            )


def load_mafs_to_db(conn, metadata, maf_pths):
    ins = metadata.tables['gdc_per_caller'].insert()
    # Insert mutation calls by batches
    ins_batch = []
    for i, mut in enumerate(read_all_mafs(maf_pths), 1):
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
    """Define the SQLite database schema."""
    Table('gdc_raw_per_caller', metadata, *create_cols(RawMutationCall))
    Table('gdc_per_caller', metadata, *create_cols(AnnotatedMutationCall))


def main(db_pth, gdc_raw_vcfs, gdc_mafs):
    logger.info(f'Build SQLite database at {db_pth}')
    metadata = MetaData()
    db_engine = create_engine('sqlite:///' + db_pth)
    define_db_schema(metadata)
    metadata.create_all(db_engine)
    conn = db_engine.connect()

    logger.info(f'Load GDC raw mutation calls from {len(gdc_raw_vcfs)} VCFs')
    load_vcfs_to_db(conn, metadata, gdc_raw_vcfs)
    logger.info(f'Load GDC annotated mutation calls from {len(gdc_mafs)} MAFs')
    load_mafs_to_db(conn, metadata, gdc_mafs)
    logger.info(f'Group calls from the multiple callers')
    conn.execute(GROUP_GDC_ANNOTATED_CALLS_SQL)

    logger.info('Complete')


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA cache_size=-8000000")
    cursor.execute("PRAGMA temp_store=MEMORY")
    cursor.execute("PRAGMA journal_mode=MEMORY")
    cursor.close()


def setup_cli(snakemake):
    # Setup logging to file
    fh = logging.FileHandler(snakemake.log[0], 'w')
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(fh)
    log_fmt = '[%(asctime)s][%(levelname)-7s] %(message)s'
    log_formatter = logging.Formatter(log_fmt, '%Y-%m-%d %H:%M:%S')
    fh.setFormatter(log_formatter)


if __name__ == '__main__':
    setup_cli(snakemake)  # noqa
    main(
        snakemake.output['db'],             # noqa
        snakemake.input['raw_vcfs'],        # noqa
        snakemake.input['annotated_mafs']   # noqa
    )
