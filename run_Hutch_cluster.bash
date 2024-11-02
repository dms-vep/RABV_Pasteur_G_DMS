#!/bin/bash

snakemake -s dms-vep-pipeline-3/Snakefile --unlock

snakemake -j 8 -s dms-vep-pipeline-3/Snakefile --software-deployment-method conda --rerun-incomplete
