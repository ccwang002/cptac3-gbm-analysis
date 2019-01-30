from collections import namedtuple
import logging
import gzip

logger = logging.getLogger(__name__)

# The class to store one sequence from FASTA
FastaEntry = namedtuple('FastaEntry', 'name seq')


def read_fasta(pth):
    """
    Reads a FASTA file.

    Each sequence is stored as a FastaEntry(name, seq).
    Note that the sequence may contain newlines.
    """
    current_entry_name = None
    current_entry_val = []
    with gzip.open(pth, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                # If this is not the first entry, we yield the previous record
                if current_entry_name is not None:
                    yield FastaEntry(current_entry_name, b''.join(current_entry_val))
                # Get the entry name by excluding the '>' symbol at the
                # beginning
                current_entry_name = line[1:].rstrip()
                # Create new empty entry
                current_entry_val = []
            else:
                # Keep appending sequence into current cDNA entry
                current_entry_val.append(line)
        yield FastaEntry(current_entry_name, b''.join(current_entry_val))


def main(cdna_fa_pth, tx_list_pth, out_fa_pth, missing_tx_pth):
    # Read the selected TX ids
    selected_tx_ids = set()
    with open(tx_list_pth, 'rb') as f:
        selected_tx_ids.update(line.rstrip(b'\n') for line in f)
    logger.info(f'Selected {len(selected_tx_ids):,d} transcript IDs from {tx_list_pth}.')

    # Subset the cDNA to a new FASTA file
    exported_cdnas = set()
    with gzip.open(out_fa_pth, 'wb', compresslevel=9) as f:
        for i, cdna in enumerate(read_fasta(cdna_fa_pth), 1):
            # The double split here makes sure it works for
            # both Ensembl and GENCODE fasta files.
            tx_id = cdna.name.split(b' ', 1)[0].split(b'|', 1)[0]

            if tx_id in selected_tx_ids:
                # Write the cDNA sequences of the selected transcripts
                # Note that we only use the transcript ID as the name
                f.write(b'>' + tx_id + b'\n' + cdna.seq)
                exported_cdnas.add(tx_id)

    logger.info(f'Finish reading {i:,d} transcript cDNAs from {cdna_fa_pth}')
    logger.info(f'{len(exported_cdnas):,d}/{len(selected_tx_ids):,d} selected transcripts '
                f'were exported to {out_fa_pth}')

    # Export the transcript IDs not found in the cDNA
    tx_ids_missing_cdna = selected_tx_ids - exported_cdnas
    with open(missing_tx_pth, 'bw') as f:
        for tx_id in sorted(tx_ids_missing_cdna):
            f.write(tx_id + b'\n')

    logger.info(f'Exported the {len(tx_ids_missing_cdna):,d} transcript IDs '
                f'that were not found in the cDNA file to {missing_tx_pth}')


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
        snakemake.input['cdna_fa'],  # noqa
        snakemake.input['tx_list'],  # noqa
        snakemake.output['out_fa'],  # noqa
        snakemake.output['missing_tx_list'],  # noqa
    )
