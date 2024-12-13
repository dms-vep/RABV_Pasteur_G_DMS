# Scratch notebooks
### Arjun Aditham, Bloom Lab (2024)

This directory contains notebooks that run scripts for custom analyses. They do not run as part of the pipeline. Some rely on data that is produced by the pipeline. Others were scripted in order to produce pre-requisite information for the pipeline. These notebooks have been edited slightly to increase readability and to update any directories that rely on additional data files.

## Notebooks that produce figures for the manuscript

Some notebooks were used to generate figures for the manuscript. Below are the notebooks and the raw figures they produced (or the analysis used to generate other figures)

'240810_LuciferaseValidations.ipynb' -- Analyzes data and produces the Figure 1D.

'AA_sub_distribution_plots_FigsS1CD.ipynb' -- produces the raw numbers for Fig S1C and the distribution plots for S1D. 

'241126_GetContacts.ipynb' -- produces list of sites within 4 Angstroms of any atom in antibody for structures of complexes are available. These sites are used to generate heatmaps of escape at specific sites.

'17C7_EscapeValidations.ipynb' -- This notebook produces Fig 3E and Fig S5A. Analyzes neutralization data for 17C7 validation assays. 

'231014_NeutCurves.ipynb' -- Produces the neut curves seen in Figure S5B.

'231027_NeutCurves.ipynb' -- Produces the neut curves seen in Figure S5B.

'241119_Heatmap_AntibodyConctact.ipynb' -- produces the heatmaps used in Figure 5. This notebook was hard-coded with sites identified using the distance cut-off rules as defined in its partner notebook 241126_GetContacts.ipynb. 

'natural_sequence_antibody_escape-PaperFigures.ipynb' -- Reformats the antibody escape versus frequency plots for Figure 6A. Produces top escape mutations listed in Figure 6C.

'241003_CircStrain-QEJ_Neuts.ipynb' -- Produces neutralization assay data for Figure 6D for the NY-2011-1202O rabies G pseudovirus. Note that we use 'QEJ' as a convenient shorthand notation for the strain.

'240825_CircStrain-AGN_Neuts.ipynb' -- Produces neutralization assay data for Figure 6D for the A12_2718 rabies G pseudovirus. Note that we use 'AGN' as a convenient shorthand notation for the strain.

'AGN-mutants_Neuts.ipynb' -- Produces Figure S6B and analyzes neutralization assay data for mutations in the paper. Note that we use 'AGN' as a convenient shorthand for the A12_2718 full strain name. 


## Notebooks

glycoprotein_positions.ipynb -- produces various csv files for pipeline, including site_numbering_map.csv, PacBio_runs.csv, and designed_mutations.csv. Note that these files must reside within the data directory for the pipeline to work as-is. 

color_pdb_functionaleffects.ipynb -- produces pdb file with b-factors updated to contain mean functional effects on cellular entry for visualizing effects on PyMOL. This notebook requires paths to a pdb file of the structure of interest and to the csv file with calculated functional effects.

update_barcode_runs.ipynb -- produces the updated barcode_runs.csv document that is input into the dms-vep-pipeline-3 workflow for analyzing functional and antibody selections. Briefly, this notebook takes in the prior barcode_runs.csv file, archives the prior set of selections into a separate user-defined document, and appends additional runs and updates the file.

