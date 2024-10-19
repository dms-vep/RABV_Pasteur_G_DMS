#!/bin/bash

# Description:
# bash script to run snakemake pipeline on cluster

# Author:
# Caleb Carr

# Parameters:
# N/A 

# stop on errors
set -e

# Signal that snakemake is running
echo "Running snakemake..."

# Run the main analysis on `slurm` cluster
snakemake \
    -j 1 \
    --cores 16 \
    --software-deployment-method conda \
    --conda-frontend conda \
    --latency-wait 60 \
    --restart-times 2 \
    --rerun-triggers mtime

# Signal that snakemake has complete
echo "Run of snakemake complete."