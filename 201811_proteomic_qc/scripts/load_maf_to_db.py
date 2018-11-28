import argparse
import logging
from pathlib import Path
from typing import NamedTuple
from sqlalchemy import (
    create_engine, event,
    MetaData, Table, Column, Integer, Text, Index,
    UniqueConstraint
)
from sqlalchemy.engine import Engine
from maf_utils import GDCMAF

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
    group_concat(caller, '|') AS callers,
    group_concat("filter", '|') AS "filter",
    variant_type, variant_classification,
    symbol, gene, transcript_id, tsl, biotype, consequence, canonical,
    hgvsc, hgvsp, hgvsp_short,
    all_effects
FROM gdc_raw_per_caller
WHERE filter IN ('PASS', 'panel_of_normals')
GROUP BY sample, chromosome, start, end, ref_allele, alt_allele
"""


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


def define_db_schema(metadata):
    columns = []
    for field, f_type in AnnotatedMutationCall._field_types.items():
        if f_type is int:
            col = Column(field, Integer())
        else:
            col = Column(field, Text())
        columns.append(col)

    Table('gdc_raw_per_caller', metadata, *columns)


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
                sample=sample, chromosome=m.chromosome, start=int(m.start_position), end=int(m.end_position),
                ref_allele=m.reference_allele, alt_allele=m.tumor_seq_allele2,
                filter=m.filter, caller=m.caller,
                t_depth=to_int(m.t_depth), t_ref_count=to_int(m.t_ref_count), t_alt_count=to_int(m.t_alt_count),
                n_depth=to_int(m.n_depth), n_ref_count=to_int(m.n_ref_count), n_alt_count=to_int(m.n_alt_count),
                variant_type=m.variant_type, variant_classification=m.variant_classification,
                symbol=m.symbol, gene=m.gene, transcript_id=m.transcript_id, tsl=m.tsl,
                biotype=m.biotype, consequence=m.consequence, canonical=m.canonical,
                hgvsc=m.hgvsc, hgvsp=m.hgvsp, hgvsp_short=m.hgvsp_short,
                all_effects=m.all_effects
            )


def load_mafs_to_db(conn, metadata, maf_pths):
    ins = metadata.tables['gdc_raw_per_caller'].insert()
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


def main(db_pth, gdc_mafs):
    logger.info(f'Load mutation calls to SQLite at {db_pth}')
    logger.info(f'... loading {len(gdc_mafs)} MAFs')
    metadata = MetaData()
    db_engine = create_engine('sqlite:///' + db_pth)

    define_db_schema(metadata)
    metadata.create_all(db_engine)

    conn = db_engine.connect()
    load_mafs_to_db(conn, metadata, gdc_mafs)

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


if __name__ == '__main__':
    setup_cli()
    main(snakemake.output['db'], snakemake.input['mafs'])
