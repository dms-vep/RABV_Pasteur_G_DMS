#!/bin/bash
#
#SBATCH -c 32

snakemake -j 32 -s dms-vep-pipeline-3/Snakefile --use-conda --rerun-incomplete