## 201811 Proteomic QC

### Setup
Also need perl 5.

    conda install snakemake-minimal


### Run
Process and load GDC mutations:

    snakemake load_gdc_mutation_to_db


### Analysis

1. `mutation_frequency.Rmd`: overall mutation rate of the samples
2. `mutation_genes.Rmd`: investigate the mutated genes and their impact
3. `idh1_idh2_hotspots.Rmd`: focus on IDH1/2 hotspot mutations
