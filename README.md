# CPTAC3 GBM analysis
The repository contains all the additional preprocessing scripts/pipelines, and
analysis scripts to generate all the results and figures for CPTAC3
Giloblastoma (GBM) cohort.

For full information about this project, please refer to [the lab wiki][wiki]
for details.

[wiki]: https://confluence.ris.wustl.edu/pages/viewpage.action?pageId=37130883



## Matt's GDC catalog
`matt_catalog` is the [GDC data catelog][matt-catelog-github] maintained by
Matt in another GitHub repository. To update the catelog to the latest version,
run git subtree at the root folder of this git repository:

```
git subtree pull --squash \
    --prefix matt_catalog \
    git@github.com:ding-lab/CPTAC3.catalog.git master
```

[matt-catelog-github]: https://github.com/ding-lab/CPTAC3.catalog



## Annotations
`annotations` is a folder storing commonly used genomic anotations by the projects.


## Consolidating discovery cohort
### `201907_locate_data` on denali



## Ad-hoc cohort
### `201904_locate_adhoc_data` on denali

### `201904_harmonize_dna_seq` on MGI
Harmonize DNA-seq (WGS and WXS) ourselves from the raw sequencing data using
the GDC workflow.

### `201904_gene_quantification` on katmai and MGI
Generate the transcript quantification (TPM) and gene-level expression of the
cohort using Salmon.

### `201904_mutation` on denali
Load the somatic mutations  of the adhoc cohort generated by WashU
(SomaticWrapper), WashU (TinDaisy) and GDC.

### `201905_somatic_wgs_cnv` on katmai
Collect the somatic CNV calls from WGS.



## Discovery cohort
### `201901_locate_discovery_data` on denali
Locate the up-to-date GDC and proteomic data release of the discovery cohort
and the local file maps on denali, katmai, and MGI.


### `201901_gene_quantification` on katmai
Generate the transcript quantification (TPM) and gene-level expression of the
discovery cohort using Salmon.


### `201901_mutation_and_variant` on denali
Load the somatic mutations and germline variants of the discovery cohort
generated by WashU and GDC.


### `201902_somatic_wgs_cnv` on katmai
Collect the somatic CNV calls from WGS.


### `201903_case_review` on katmai
Review some cases that may not be glioblastoma.  It also tests the somatic CNV
calling using WGS/WXS sequencing coverage, and GATK4 PoN pipeline.


### `201903_methylation_array` on denali
Annotate the processed EPIC methylation microarray.



## Proteomic QC
### `201811_locate_gdc_data` on denali
📍 *Use `201901_locate_discovery_data` instead*<br>
Locate the current GDC data release of CPTAC3 GBM cohort and their local
whereabouts.


### `201811_proteomic_qc` on denali
The preliminary analysis to help Tao from PNNL to design the proteomic
experiment plexes. The main goal is to address the 5 tumors that are clustered
with GTEx normals using single-shot mass spectrometry.


### `201811_rnaseq_salmon` on katmai
📍 *Use `201901_gene_quantification` instead*<br>
RNA quantification of high-confident protein coding trasncripts (Ensembl v90).
Only transcripts of protein coding with TSL 1 or TSL NA (e.g. single exon
trasncripts) were included. Since the mapping rate using this approach is low,
we also explore some other RNA-seq pipelines:

- Salmon using all transcripts
- Alignment-based pipeline: STAR + htseq


### `201811_wxs_bam_readcount` on katmai
WXS BAM readcount on regions of interest.
