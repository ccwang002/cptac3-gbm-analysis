rsync -az --info=progress2 --exclude-from=scripts/cnv_rsync_exclude_patterns.list \
    katmai:/diskmnt/Projects/cptac_downloads_4/CPTAC3_GBM/201905_somatic_wgs_cnv/results.GBM-not-in-discover.katmai/ \
    external_data/bicseq2_cnv

rsync -az --info=progress2 --exclude-from=scripts/cnv_rsync_exclude_patterns.list \
    katmai:/diskmnt/Projects/cptac_downloads_4/CPTAC3_GBM/201902_somatic_wgs_cnv/run_cases.GBM-subset/ \
    external_data/bicseq2_cnv

rsync -az --info=progress2 --exclude-from=scripts/cnv_rsync_exclude_patterns.list \
    katmai:/diskmnt/Projects/cptac_downloads_4/CPTAC3_GBM/201908_somatic_wgs_cnv/results.GBM-GDC-BAM-20190802.katmai/ \
    external_data/bicseq2_cnv
