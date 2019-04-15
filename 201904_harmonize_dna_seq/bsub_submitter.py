#!/usr/bin/env python3
# From https://github.com/broadinstitute/viral-ngs/blob/master/pipes/Broad_LSF/run-pipe.sh
import os
import sys
import re
from snakemake.utils import read_job_properties

lsf_log_dir = sys.argv[-2]
jobscript = sys.argv[-1]
props = read_job_properties(jobscript)

# set up job name, project name
jobname = "{rule}_job{jobid}".format(rule=props["rule"], jobid=props["jobid"])

# -E is a pre-exec command, that reschedules the job if the command fails
#   in this case, if the data dir is unavailable (as may be the case for a hot-mounted file path)
# cmdline = f'bsub -J {jobname} -r -E "ls {DATADIR}" '
cmdline = f'bsub -J {jobname} -r '

cluster = props['cluster']
if 'docker_image' in cluster:
    docker_image = cluster['docker_image']
    cmdline += f'-a "docker({docker_image})" '

# Add queue name
if 'queue' in cluster:
    queue = cluster['queue']
    cmdline += f'-q {queue} '

# Log output to a file
if "-N" not in props["params"].get("LSF", ""):
    cmdline += f"-oo {lsf_log_dir}/LSF_{jobname}.log "

# Pass memory and cpu resource request to LSF
ncpus = props['threads']
mem = props.get('resources', {}).get('mem')
if mem:
    cmdline += (
        f'-R "select[maxmem>{mem} && ncpus>={ncpus}] rusage[mem={mem}]" '
        f'-n {ncpus} -M {mem}000 '
    )
else:
    cmdline += (
        f'-R "select[ncpus>={ncpus}]" '
        f'-n {ncpus} '
    )

# Add rule-specific LSF parameters (e.g. queue, runtime)
cmdline += props["params"].get("LSF", "") + " "

# figure out job dependencies
dependencies = sys.argv[1:-2]
if dependencies:
    # Create the LSF dependency expression
    dep_expr = " && ".join(dependencies)
    cmdline += f"-w '{dep_expr}' "

# the actual job
cmdline += jobscript

# the part that strips bsub's output to just the job id
cmdline += r" | tail -1 | cut -f 2 -d \< | cut -f 1 -d \>"

# call the command
os.system(cmdline)

