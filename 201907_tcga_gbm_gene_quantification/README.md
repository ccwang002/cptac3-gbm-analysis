# 201907_tcga_gbm_gene_quantification

The files under `tcga_gbm_gdc_manifest` were generated on GDC Portal using the following query:

    cases.project.project_id in ["TCGA-GBM"]
    and files.analysis.metadata.read_groups.read_length >= 76
    and files.analysis.metadata.read_groups.read_length <= 76
    and files.data_category in ["Sequencing Reads"]
    and files.data_format in ["BAM"]
    and files.experimental_strategy in ["RNA-Seq"]
