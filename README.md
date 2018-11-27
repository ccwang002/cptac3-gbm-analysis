## CPTAC3 GBM analysis
The repository contains all the additional preprocessing scripts/pipelines, and
analysis scripts to generate all the results and figures for CPTAC3
Giloblastoma (GBM) cohort.

For full information about this project, please refer to [the lab wiki for details][wiki].

[wiki]: https://confluence.ris.wustl.edu/pages/viewpage.action?pageId=37130883


### `201811_locate_gdc_data`
Locate the current GDC data release of CPTAC3 GBM cohort and their local whereabouts.


### `201811_proteomic_qc`
The preliminary analysis to help Tao from PNNL to design the proteomic
experiment plexes. The main goal is to address the 5 tumors that are clustered
with GTEx normals using single-shot mass spectrometry.


### `201811_rnaseq_salmon`
RNA quantification of high-confident protein coding trasncripts (Ensembl v90). Only transcripts of protein coding with TSL 1 or TSL NA (e.g. single exon trasncripts) were included. Note that CPTAC used total RNA-seq, so the mapping rate using this approach is low.
