# 201901_locate_discovery_data
This project finds the hg38 genomic data from GDC, including:

- WXS BAMs
- WGS BAMs
- RNA-seq FASTQs (R1 and R2)
- miRNA-seq FASTQ


### Notebooks
1. `locate_discovery_samples.Rmd`: Find the GDC UUIDs of the discovery samples.
2. `genome_data_local_availability.Rmd`: Check the genomic data availability on the lab servers.
3. `generate_gdc_vcf_manifest.Rmd`: Generate the GDC manifest of raw and annotated VCFs of the discovery samples.


### Annotations
Folder `annotations` stores various annotation (R) objects that will be shared for various pipeline.

- `seqinfo_GRCh38.d1.vd1.rds`: GDC genome information (`GRCh38.d1.vd1`) as a SeqInfo R object.
- `cytoband_hg38.rds`: cytoband information (canonical hg38 chromosomes only) as a tibble R object.
- `EnsDb.Hsapiens.v94.sqlite`: Ensembl v94 annotation as an ensembldb R object downloaded from [AnnotationHub][ensdb].
- `ucscToEnsembl.txt.gz`: Chromosome name conversion between UCSC and Ensembl; downloaded from [UCSC][ucsc-chrom].

[ensdb]: http://s3.amazonaws.com/annotationhub/AHEnsDbs/v94/EnsDb.Hsapiens.v94.sqlite
[ucsc-chrom]: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/ucscToEnsembl.txt.gz
