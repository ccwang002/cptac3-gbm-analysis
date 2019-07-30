# 201907_locate_discovery_data
This project tries to locate the data of GBM discovery cohort.

## Notebooks
1. `notebooks/locate_all_gbm_samples.Rmd`: Locate all the available GBM cases with sequencing data
2. `notebooks/download_raw_seq_data.Rmd`: Download newer/changed sequencing data from GDC
3. `notebooks/find_cases_to_run_mutation.Rmd`: Find new cases to run somatic mutation calling
4. `notebooks/download_mirna_fq.Rmd`: Download more miRNA-seq FASTQs

Don't run (as now the data have been updated):
- `notebooks/gen_unharmonized_case_list.Rmd`: Generate unharmonized cases (in GDC harmonization batch 3)