# 201907_gdc_rnaseq_counts
Process the GDC RNA-Seq counts data for TCGA GBM and CPTAC GBM.

TCGA GBM counts files were selected based on the following GDC query on GDC Portal:

    cases.project.project_id in ["TCGA-GBM"] 
    and files.data_category in ["Transcriptome Profiling"] 
    and files.experimental_strategy in ["RNA-Seq"]

CPTAC GBM counts files were filtered based on Mathangi's manifest (not released yet)