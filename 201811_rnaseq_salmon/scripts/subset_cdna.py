from collections import namedtuple
import gzip

FastaEntry = namedtuple('FastaEntry', 'name seq')


def read_fasta(pth):
    """Reads a FASTA file"""
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

# Read transcript IDs
selected_tx_ids = []
with open(snakemake.input['tx_list']) as f:
    selected_tx_ids = [tx.strip().encode() for tx in f]


with gzip.open(snakemake.output[0], 'wb', compresslevel=5) as f:
    for cdna in read_fasta(snakemake.input['cdna_fa']):
        tx_id = cdna.name.split(b' ', 1)[0]
        if tx_id in selected_tx_ids:
            f.write(b'>' + cdna.name + b'\n' + cdna.seq)
            f.flush()
