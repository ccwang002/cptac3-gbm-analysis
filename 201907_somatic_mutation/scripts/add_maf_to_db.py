import gzip
import logging
from typing import NamedTuple
from sqlalchemy import (
    create_engine, event,
    MetaData, Table, Column, Integer, Text,
)
from sqlalchemy.engine import Engine
from maf_utils import TrailingTabTrimmedMAF

logger = logging.getLogger(__name__)
BATCH_SIZE = 5000


class MAFMutationCall(NamedTuple):
    sample: str
    chromosome: str
    start: int
    end: int
    ref_allele: str
    alt_allele: str
    filter: str
    callers: str
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
    variant_class: str
    consequence: str
    canonical: str
    hgvsc: str
    hgvsp: str
    hgvsp_short: str
    all_effects: str

    @classmethod
    def create_cols(cls):
        """Create SQLAlchemy columns based on a named tuple class."""
        columns = []
        for field, f_type in cls._field_types.items():
            if f_type is int:
                col = Column(field, Integer())
            else:
                col = Column(field, Text())
            columns.append(col)
        return columns


    @classmethod
    def define_db_schema(cls, metadata, db_table_name):
        """Define the SQLite database schema."""
        Table(db_table_name, metadata, *cls.create_cols())


def read_maf(maf_pth, callers_pth):
    """Read a MAF and the associated callers file."""
    def to_int(x):
        if x == '':
            return 0
        elif x == '.':
            return None
        else:
            return int(x)

    maf = TrailingTabTrimmedMAF(maf_pth)
    callers_per_mut = gzip.open(callers_pth, 'rt')
    assert next(callers_per_mut).rstrip('\n') == 'callers'

    for m, callers in zip(maf, callers_per_mut):
        sample = m.tumor_sample_barcode[:-2]
        yield MAFMutationCall(
            sample=sample, chromosome=m.chromosome,
            start=int(m.start_position), end=int(m.end_position),
            ref_allele=m.reference_allele, alt_allele=m.tumor_seq_allele2,
            filter=m.filter, callers=callers.rstrip('\n'),
            t_depth=to_int(m.t_depth), t_ref_count=to_int(m.t_ref_count),
            t_alt_count=to_int(m.t_alt_count),
            n_depth=to_int(m.n_depth), n_ref_count=to_int(m.n_ref_count),
            n_alt_count=to_int(m.n_alt_count),
            variant_type=m.variant_type, variant_classification=m.variant_classification,
            symbol=m.symbol, gene=m.gene, transcript_id=m.transcript_id, tsl=m.tsl,
            biotype=m.biotype, variant_class=m.variant_class, consequence=m.consequence,
            canonical=m.canonical,
            hgvsc=m.hgvsc, hgvsp=m.hgvsp, hgvsp_short=m.hgvsp_short,
            all_effects=m.all_effects
        )


def load_maf_to_db(conn, metadata, maf_pth, callers_pth, db_table_name):
    """Load MAF to a given table."""
    ins = metadata.tables[db_table_name].insert()
    # Insert mutation calls by batches
    ins_batch = []
    for i, mut in enumerate(read_maf(maf_pth, callers_pth), 1):
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



def main(db_pth, maf_pth, callers_pth, db_table_name):
    logger.info(f'Add new table {db_table_name} to the SQLite database at {db_pth}')
    metadata = MetaData()
    db_engine = create_engine(f'sqlite:///{db_pth}')
    MAFMutationCall.define_db_schema(metadata, db_table_name)
    metadata.create_all(db_engine)
    conn = db_engine.connect()

    logger.info(f'Load MAF from {maf_pth} to table {db_table_name}')
    load_maf_to_db(conn, metadata, maf_pth, callers_pth, db_table_name)

    logger.info('Complete')


if __name__ == '__main__':
    setup_cli(snakemake)  # noqa
    main(
        snakemake.params['db'],   # noqa
        snakemake.input['maf'],   # noqa
        snakemake.input['callers'],         # noqa
        snakemake.params['db_table_name'],  # noqa
    )
