# Data

This directory contains input data for multiple steps in the `dms-vep-3` pipeline. These files are detailed below, sorted by analytical step.

## PacBio sequencing analysis to link barcodes with glycoprotein variants

[PacBio_amplicon.gb](PacBio_amplicon.gb): Genbank file of the PacBio amplicons with annotated features used in [alignparse](https://jbloomlab.github.io/alignparse/). Must have *gene* (the gene of interest) and *barcode* features. Regions of the glycoprotein ORF that are not mutagenized are separately annotated, for example with *gene_flank5*. 

[PacBio_feature_parse_specs.yaml](PacBio_feature_parse_specs.yaml): How to parse the PacBio amplicon using [alignparse](https://jbloomlab.github.io/alignparse/). Note this has been edited to contain an extra feature *gene_flank5* and *gene _flank3*. These region of the glycoprotein were not mutagenized, so any mutations here should be ignored in analysis.

[PacBio_runs.csv](PacBio_runs.csv): List of PacBio CCS FASTQs used to link barcodes to variants. This file contains the following columns
 - `library`: name of the library sequenced, A or B (replicate libraries).
 - `run`: Library name with run date appended. For example, 
    + *A_240124*: Library A PacBio sequencing submitted on 24-Jan-2024.
    + *B_240124*: Library B PacBio sequencing submitted on 24-Jan-2024.
    + *B_240212*: Library B PacBio re-sequencing submitted on 12-Feb-2024.
 - `fastq`: FASTQ file from running CCS
    + We asked for genomics to convert bam files into fastq on our behalf, though `samtools` is an alternative way to convert bam files to fastq for analysis.
  
[neutralization_standard_barcodes.csv](neutralization_standard_barcodes.csv): These are a minimal set of barcodes *reserved* as sequencing standards in antibody selections. These are barcodes in a VSV-G pseudotyped lentiviral particles containing barcoded genomes encoding mCherry.  
  
[site_numbering_map.csv](site_numbering_map.csv): Maps sequential (1, 2, ...n) numbering of positions to a "reference" numbering scheme for analysis. This file contains the following columns:`
 - `sequential_site`: site within mutagenized region of glycoprotein, numbered from 1 onwards. While `reference_site` does not have to begin with 1, `sequential_site` must start at 1.
 - `sequential_wt`: wild type amino acid residue at this position
 - `reference_site`: numbering used (conventional in field) or within entire ORF in case only part of glycoprotein is subject to DMS.
 - `reference_wt`: wild type amino acid residue within reference strain
 - `region`: domain of glycoprotein where position residues, eg. *signal peptide*, *ectodomain*. 

[designed_mutations.csv](designed_mutations.csv) lists the designed amino-acid mutations and stop-codon mutations, using sequential numbering of the mutagenized region.

[./gene_sequence/](gene_sequence): This directory contains files which contain the protein and nucleic acid sequence of the glycoprotein. Note that these files only include the mutagenized ectodomain region of the glycoprotein. 

## Illumina sequencing analysis of functional and antibody selections

[barcode_runs.csv](barcode_runs.csv) contains the following columns:
 - `sample`: sample name, must be unique among barcode runs. Sample name must begin with `<library>-<YYMMDD>` where `<library>` is the library (A or B, for example) and `<YYMMDD>` is the date of the selection. A recommended format is `<library>-<YYMMDD>-<description>-<replicate>`. `<description>` is a string description with underscores but no dashes (eg. *VSVG* is for VSV-G virus selections for barcode normalization, *no_antibody* means functional selection only). `<replicate>` is a number.
 - `library`: name of library, must match a library in the barcode-variant table
 - `date`: date of experimental run, specified in a format parseable to a date by `pandas`.
 - `fastq_R1`: path to one more FASTQ R1 sequencing files. Multiple sequencing files should be delimited using a ';' character in between file paths.

[func_effects_config.yml](func_effects_config.yml) has the configuration for analyzing functional effects of mutations.

[antibody_escape_config.yml](antibody_escape_config.yml) has the configuration for analyzing effects of mutations on escape from antibodies or sera.


## Visualization on known structures of glycoprotein

The below structures were used to project functional effects onto structures for mechanistic rationalization. This is now in the [./PDB_Structures/](PDB_Structures) subdirectory.

7u9g.pdb: Cryo-EM structure of trimeric, pre-fusion Rabies glycoprotein (Pasteur strain) published in [Callaway, et al, *Sci Adv* 2022](https://www.science.org/doi/10.1126/sciadv.abp9151). This structure contains the Fab fragment of the antibody RVA122 bound. 

6lgx.pdb: Structure of extended intermediate state of Rabies glycoprotein during fusion process, as published in [Yang, et al, *Cell Host Microbe*, 2020](https://www.sciencedirect.com/science/article/pii/S1931312819306419?via%3Dihub).

6lgw.pdb: basic pH structure of rabies glycoprotein, companion structure with 6lgx.pdb also published in [Yang, et al, *Cell Host Microbe*, 2020](https://www.sciencedirect.com/science/article/pii/S1931312819306419?via%3Dihub).

8a1e.pdb: Structure of 17C7 antibody bound to rabies glycoprotein, as published in [Ng, et al., *Cell Host Microbe*, 2022](https://doi.org/10.1016/j.chom.2022.07.014).

6tou.pdb: Structure of pH domain of rabies glycoprotein (fragment) bound to RVC20, as published in [Hellert, et al., *Nat Comms*, 2020](https://doi.org/10.1038/s41467-020-14398-7).

A configuration YAML that helps with viewing these structures with `dms-viz` is in [dms_viz_config.yml](dms_viz_config.yml).
