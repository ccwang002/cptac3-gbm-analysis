# Annotations
This folder stores various annotation objects that are re-used by many projects.

- `seqinfo_GRCh38.d1.vd1.rds`: GDC genome information (`GRCh38.d1.vd1`) as a SeqInfo R object.
- `cytoband_hg38.rds`: cytoband information (canonical hg38 chromosomes only) as a tibble R object.
- `EnsDb.Hsapiens.v94.sqlite`: Ensembl v94 annotation as an ensembldb R object downloaded from [AnnotationHub][ensdb-v94].
- `EnsDb.Hsapiens.v96.sqlite`: Ensembl v96 annotation as an ensembldb R object downloaded from [AnnotationHub][ensdb-v96].
- `ucscToEnsembl.txt.gz`: Chromosome name conversion between UCSC and Ensembl; downloaded from [UCSC][ucsc-chrom].
- `hg38_exome_targeted_regions.bed.gz`:
    CPTAC3 whole exome sequencing (WXS) target regions in hg38.
    Originally from `~/Box/Ding_Lab/Projects_Current/CPTAC3-GBM/Resources/Annotations/WXS/hg38_exome_targeted_regions.bed.gz`.
- `Batch1through6_sample_attributes.csv.gz`: CPTAC genomic sample to case mapping table got from Mathangi.

[ensdb-v94]: http://s3.amazonaws.com/annotationhub/AHEnsDbs/v94/EnsDb.Hsapiens.v94.sqlite
[ensdb-v96]: http://s3.amazonaws.com/annotationhub/AHEnsDbs/v96/EnsDb.Hsapiens.v96.sqlite
[ucsc-chrom]: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/ucscToEnsembl.txt.gz
