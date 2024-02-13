# Scratch notebooks
### Arjun Aditham, Bloom Lab (2024)

This directory contains notebooks that run scripts for custom analyses. They do not run as part of the pipeline. Some rely on data that is produced by the pipeline. Others were scripted in order to produce pre-requisite information for the pipeline.

## Notebooks

glycoprotein_positions.ipynb -- produces various csv files for pipeline, including site_numbering_map.csv, PacBio_runs.csv, and designed_mutations.csv. Note that these files must reside within the data directory for the pipeline to work as-is. 


color_pdb_functionaleffects.ipynb -- produces pdb file with b-factors updated to contain mean functional effects on cellular entry for visualizing effects on PyMOL. This notebook requires paths to a pdb file of the structure of interest and to the csv file with calculated functional effects.

update_barcode_runs.ipynb -- produces the updated barcode_runs.csv document that is input into the dms-vep-pipeline-3 workflow for analyzing functional and antibody selections. Briefly, this notebook takes in the prior barcode_runs.csv file, archives the prior set of selections into a separate user-defined document, and appends additional runs and updates the file.