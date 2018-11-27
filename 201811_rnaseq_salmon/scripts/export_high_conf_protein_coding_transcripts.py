from gffutils import FeatureDB

db = FeatureDB(snakemake.input[0])

tsl_good_txes = []

for tx in db.all_features(featuretype='transcript'):
    # Keep only protein coding transcripts
    if tx.attributes['gene_biotype'][0] != 'protein_coding':
        continue
    tx_ver = tx.attributes['transcript_version'][0]
    try:
        tsl = tx.attributes['transcript_support_level'][0]
        if tsl == '1' or tsl == 'NA':
            tsl_good_txes.append((tx.id, tx_ver))
    except KeyError:
        # The TSL is considered NA as well
        tsl_good_txes.append((tx.id, tx_ver))


with open(snakemake.output[0], 'w') as f:
    for tx_id, tx_ver in tsl_good_txes:
        f.write(f'{tx_id}.{tx_ver}\n')

