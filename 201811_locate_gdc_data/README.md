## 201811_locate_gdc_data
This project finds the hg38 genomic data from GDC, including:

- WXS BAMs
- WGS BAMs
- RNA-seq FASTQs (R1 and R2)
- miRNA-seq FASTQ

### `notebooks`
The R notebooks produces various sample list, subsetted manifests, and sample-wise GDC UUIDs.

1. `locate_all_GBM_samples.Rmd`: Generate the genomic sample list of 60 cases; Summary the GDC UUIDs.
2. `generate_manifest_for_proteomic_qc.Rmd`: Generate the manifest files for GDC raw and annotated VCFs.


### `matt_catelog`
`matt_catelog` is the [GDC data catelog][matt-catelog-github] maintained by
Matt available on GitHub. To update the catelog to the latest version, run git
subtree at the root:

```
git subtree pull --squash \
    --prefix 201811_locate_gdc_data/matt_catelog \
    git@github.com:ding-lab/CPTAC3.catalog.git master
```

[matt-catelog-github]: https://github.com/ding-lab/CPTAC3.catalog
