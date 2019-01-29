# 201901_locate_discovery_data
This project finds the hg38 genomic data from GDC, including:

- WXS BAMs
- WGS BAMs
- RNA-seq FASTQs (R1 and R2)
- miRNA-seq FASTQ


### `notebooks`
1. `locate_discovery_samples.Rmd`: Find the GDC UUIDs of the discovery samples.
2. `genome_data_local_availability.Rmd`: Check the genomic data availability on the lab servers.