bsub -q research-hpc -R "rusage[internet_download_mbps=100]" \
    -N -u "liang-bo.wang@wustl.edu" \
    -oo gdc_client_bsub.log \
    -a 'docker(quay.io/biocontainers/gdc-client:1.3.0--py27h1341992_3)' \
    /usr/local/bin/gdc-client download \
        -n 5 --retry-amount 2 \
        -t token/gdc-user-token.2019-07-05T17_31_20.190Z.txt \
        -m tcga_gbm_gdc_manifest/gdc_manifest.2019-07-05.txt \
        -d external_data/tcga_gbm_gdc
