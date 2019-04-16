#!/bin/bash
set -euo pipefail  # strict Bash mode
echoerr() { printf "%s\n" "$*" >&2; }

# Load the conda environment
if [ -z ${CONDA_DEFAULT_ENV+x} ]; then
    echoerr "Need to activate the conda env before running snakemake. Aborted."
    exit 1
fi

readonly master_log_pth='logs/lsf_master.log'
# Clean up the existing master log file
if [ -f $master_log_pth ]; then
    rm $master_log_pth
fi
echoerr "Launch the LSF/bsub master job to run snakemake ..."
echoerr "... master log location: $master_log_pth"
echoerr "... additional parameters: $@"
# Run snakemake via bsub
bsub -a "docker(lbwang/dailybox)" -q research-hpc \
    -N -u 'liang-bo.wang@wustl.edu' \
    -o  $master_log_pth \
    env SHELL="/bin/bash" \
    snakemake \
        --rerun-incomplete --nolock \
        --latency-wait 15 \
        --cluster "$PWD/bsub_submitter.py {dependencies} lsf_logs" \
        --cluster-config bsub_config.json \
        -p --reason --verbose \
        "$@"
