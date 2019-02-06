import logging
from typing import NamedTuple
from sqlalchemy import (
    create_engine,
    MetaData, Table
)
from maf_utils import WashUMAF
from add_gdc_maf import create_cols, setup_cli, BATCH_SIZE

logger = logging.getLogger(__name__)


class WashUMutationCall(NamedTuple):
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
    consequence: str
    canonical: str
    hgvsc: str
    hgvsp: str
    hgvsp_short: str
    all_effects: str


def read_washu_maf(maf_pth):
    """Read the WashU MAF file."""
    def to_int(x):
        if x == '':
            return 0
        else:
            return int(x)

    for m in WashUMAF(maf_pth):
        sample = m.tumor_sample_barcode[:-2]
        yield WashUMutationCall(
            sample=sample, chromosome=m.chromosome,
            start=int(m.start_position), end=int(m.end_position),
            ref_allele=m.reference_allele, alt_allele=m.tumor_seq_allele2,
            filter=m.filter, callers=m.callers,
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


def define_db_schema(metadata):
    """Define the SQLite database schema."""
    Table('washu', metadata, *create_cols(WashUMutationCall))


def load_washu_maf_to_db(conn, metadata, maf_pth):
    """Load WashU MAF to the database."""
    ins = metadata.tables['washu'].insert()
    # Insert mutation calls by batches
    ins_batch = []
    for i, mut in enumerate(read_washu_maf(maf_pth), 1):
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


def main(db_pth, washu_maf_pth, samples):
    logger.info(f'Build SQLite database at {db_pth}')
    metadata = MetaData()
    db_engine = create_engine('sqlite:///' + db_pth)
    define_db_schema(metadata)
    metadata.create_all(db_engine)
    conn = db_engine.connect()

    # Only load the mutation calls of the specified samples
    samples = set(samples)
    logger.info(f'Subset the MAF to include only {len(samples)} samples')

    logger.info(f'Load WashU annotated mutation calls from {washu_maf_pth}')
    load_washu_maf_to_db(conn, metadata, washu_maf_pth)

    logger.info('Complete')


if __name__ == '__main__':
    setup_cli(snakemake)  # noqa
    main(
        snakemake.params['db'],         # noqa
        snakemake.input['washu_maf'],   # noqa
        snakemake.params['samples']     # noqa
    )
