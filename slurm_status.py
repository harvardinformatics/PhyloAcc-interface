#!/usr/bin/python
#############################################################################
# Script to use with snakemake to check slurm status of jobs
# Modified from: https://github.com/snakemake/snakemake/issues/759
# with elements from: https://github.com/Snakemake-Profiles/slurm/blob/master/%7B%7Bcookiecutter.profile_name%7D%7D/slurm-status.py
#
# Gregg Thomas, August 2021
#############################################################################

import subprocess
import sys

#############################################################################

jobid = sys.argv[-1]
# SLURM job id as input parameter

running_status = ["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"];
# Statuses that indicate the job is running
# "COMPLETED" indicates the job is finished
# All other statuses indicate the job has failed

STATUS_ATTEMPTS = 20;
# Number of times to try to check the status 

output = "";

for i in range(STATUS_ATTEMPTS):
    try:
        output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip());
        break;
    except:
        continue;
    # Check job status with sacct
# Try to check status multiple times

if output != "":
    if "COMPLETED" in output:
        print("success");
    elif any(r in output for r in running_status):
        print("running");
    else:
        print("failed");
# Compare the output of sacct with the statuses provided as "running" statuses
