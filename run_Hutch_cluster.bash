#!/bin/bash

snakemake -j 16 -s dms-vep-pipeline-3/Snakefile --software-deployment-method conda --rerun-incomplete
