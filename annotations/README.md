# Annotations
This folder stores various annotation objects that are re-used by many projects.

Overall genome information:
- `seqinfo_GRCh38.d1.vd1.rds`: GDC genome information (`GRCh38.d1.vd1`) as a SeqInfo R object.
- `cytoband_hg38.rds`: cytoband information (canonical hg38 chromosomes only) as a tibble R object.
- `ucscToEnsembl.txt.gz`: Chromosome name conversion between UCSC and Ensembl; downloaded from [UCSC][ucsc-chrom].

Gene annotation:
- `EnsDb.Hsapiens.v94.sqlite`: Ensembl v94 annotation as an ensembldb R object downloaded from [AnnotationHub][ensdb-v94].
- `EnsDb.Hsapiens.v96.sqlite`: Ensembl v96 annotation as an ensembldb R object downloaded from [AnnotationHub][ensdb-v96].

Specific to CPTAC:
- `hg38_exome_targeted_regions.bed.gz`: CPTAC3 whole exome sequencing (WXS) target regions in hg38. Downloaded from [Box][hg38-exome-roi].
- `Batches1through9_samples_attribute.with_disease_code.tsv.gz`: CPTAC genomic sample to case mapping table got from Mathangi; Downloaded from [Box][batch-sample-attr].

From GDC:
- `gencode.gene.info.v22.tsv`: [GDC][gdc-gene-info-tsv]'s gene annotation for RNA-seq related pipelines.


[ensdb-v94]: http://s3.amazonaws.com/annotationhub/AHEnsDbs/v94/EnsDb.Hsapiens.v94.sqlite
[ensdb-v96]: http://s3.amazonaws.com/annotationhub/AHEnsDbs/v96/EnsDb.Hsapiens.v96.sqlite
[ucsc-chrom]: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/ucscToEnsembl.txt.gz
[hg38-exome-roi]: https://wustl.box.com/s/3sxgoyrw3eoca3x8au4w1ec7n4on7alp
[batch-sample-attr]: https://wustl.box.com/s/gew721788cftwglb7xqw585erbqaglh7
[gdc-gene-info-tsv]: https://api.gdc.cancer.gov/data/b011ee3e-14d8-4a97-aed4-e0b10f6bbe82
