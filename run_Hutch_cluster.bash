#!/bin/bash

snakemake -j 32 -s dms-vep-pipeline-3/Snakefile --software-deployment-method conda --conda-frontend conda --rerun-incomplete