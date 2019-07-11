rsync -az --info=progress2 --exclude-from=scripts/cnv_rsync_exclude_patterns.list \
    vw3.gsc.wustl.edu:/gscmnt/gc2541/cptac3_analysis/gbm_somatic_wgs_cnv/results.GBM-custom-harm.MGI/ \
    external_data/bicseq2_cnv

rsync -az --info=progress2 --exclude-from=scripts/cnv_rsync_exclude_patterns.list \
    vw3.gsc.wustl.edu:/gscmnt/gc2541/cptac3_analysis/gbm_somatic_wgs_cnv/results.GBM-custom-harm-20190710.MGI/ \
    external_data/bicseq2_cnv

rsync -az --info=progress2 --exclude-from=scripts/cnv_rsync_exclude_patterns.list \
    katmai:/diskmnt/Projects/cptac_downloads_4/CPTAC3_GBM/201905_somatic_wgs_cnv/results.GBM-not-in-discover.katmai/ \
    external_data/bicseq2_cnv

rsync -az --info=progress2 --exclude-from=scripts/cnv_rsync_exclude_patterns.list \
    katmai:/diskmnt/Projects/cptac_downloads_4/CPTAC3_GBM/201902_somatic_wgs_cnv/run_cases.GBM-subset/ \
    external_data/bicseq2_cnv
