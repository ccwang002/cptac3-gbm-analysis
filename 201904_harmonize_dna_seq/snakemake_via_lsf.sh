#!/bin/bash
set -euo pipefail  # strict Bash mode
echoerr() { printf "%s\n" "$*" >&2; }

# Load the conda environment
if [ -z ${CONDA_DEFAULT_ENV+x} ]; then
    echoerr "Need to activate the conda env before running snakemake. Aborted."
    exit 1
fi

# Run snakemake via bsub
readonly master_log_pth='logs/lsf_master.log'
echoerr "Launch the LSF/bsub master job to run snakemake ..."
echoerr "... master log location: $master_log_pth"
echoerr "... additional parameters: $@"
bsub -a "docker(lbwang/dailybox)" -q research-hpc \
    -N -u 'liang-bo.wang@wustl.edu' \
    -o  $master_log_pth \
    env SHELL="/bin/bash" \
    snakemake \
        --rerun-incomplete --nolock \
        --cluster "bsub_submitter.py {dependencies} lsf_logs" \
        --cluster-config bsub_config.json \
        -p \
        "$@"
