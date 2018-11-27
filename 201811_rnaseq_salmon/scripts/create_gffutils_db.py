import gffutils

gffutils.create_db(
    snakemake.input[0],
    dbfn=snakemake.output[0],
    verbose=True,
    merge_strategy='error',
    disable_infer_transcripts=True,
    disable_infer_genes=True,
)
