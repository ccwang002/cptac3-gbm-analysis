# 201907_gdc_rnaseq_counts
Process the GDC RNA-Seq counts data for TCGA GBM and CPTAC GBM.

TCGA GBM counts files were selected based on the following GDC query on GDC Portal:

    cases.project.project_id in ["TCGA-GBM"]
    and files.data_category in ["Transcriptome Profiling"]
    and files.experimental_strategy in ["RNA-Seq"]

CPTAC GBM counts files were filtered based on Mathangi's manifest of GDC harmonization batch 1-3 (not released yet).


## Notebooks
- `notebooks/collect_tcga_gbm_rnaseq.Rmd`: Summarize all TCGA GBM count matrices (FPKM, FPKM-UQ, HTSeq counts)
- `notebooks/collect_cptac_gbm_rnaseq.Rmd`: Summarize all CPTAC GBM count matrices and STAR outputs
- `notebooks/parse_tcga_clinical_data.Rmd`: Parse TCGA clinical data
