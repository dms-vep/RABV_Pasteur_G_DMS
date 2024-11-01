#!/bin/bash
#
#SBATCH -c 8

snakemake -j 8 -s dms-vep-pipeline-3/Snakefile --software-deployment-method conda --rerun-incomplete
