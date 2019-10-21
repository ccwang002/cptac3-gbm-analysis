# 201910_novel_splicing
Collect novel splicing transcript from the CPTAC_splicing pipeline output


    mkdir -p external_data/summary_dats
    mkdir -p processed_data/splicing_gtfs
    parallel 'rsync -a {} external_data/summary_dats/' :::: splicing_analysis_summary.lists

Run the notebooks/collect_splicing_gtfs.Rmd. Then copy the splicing GTFs by

    cat processed_data/gbm_all_gtfs.tsv \
        | parallel -j8 --bar --colsep '\t' --header : \
            'gzip -c -9 {gtf_path} > processed_data/splicing_gtfs/{case_id}.gtf.gz'
