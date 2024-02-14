# Data

This directory contains input data for multiple steps in the 'dms-vep-3' pipeline. These files are detailed below, sorted by analytical step.

## PacBio sequencing analysis to link barcodes with glycoprotein variants

[PacBio_amplicon.gb](PacBio_amplicon.gb): Genbank file of the PacBio amplicons with annotated features used in [alignparse](https://jbloomlab.github.io/alignparse/). Must have *gene* (the gene of interest) and *barcode* features. Regions of the glycoprotein ORF that are not mutagenized are separately annotated, for example with *gene_flank5*. 

[PacBio_feature_parse_specs.yaml](PacBio_feature_parse_specs.yaml): How to parse the PacBio amplicon using [alignparse](https://jbloomlab.github.io/alignparse/). Note this has been edited to contain an extra feature *gene_flank5*. This region of the glycoprotein was not mutagenized, so any mutations here should be ignored in analysis.

[PacBio_runs.csv](PacBio_runs.csv): List of PacBio CCS FASTQs used to link barcodes to variants. This file contains the following columns
 - `library`: name of the library sequenced, A or B (replicate libraries).
 - `run`: Library name with run date appended. For example, 
    + *A_240124*: Library A PacBio sequencing submitted on 24-Jan-2024.
    + *B_240124*: Library B PacBio sequencing submitted on 24-Jan-2024.
 - `fastq`: FASTQ file from running CCS
    + We asked for genomics to convert bam files into fastq on our behalf, though 'samtools' is an alternative way to convert bam files to fastq for analysis.
  
[neutralization_standard_barcodes.csv](neutralization_standard_barcodes.csv): These are a minimal set of barcodes *reserved* as sequencing standards in antibody selections. These are barcodes in a VSV-G pseudotyped lentiviral particles containing barcoded genomes encoding mCherry.  
  
## Site numbering
[site_numbering_map.csv](site_numbering_map.csv): Maps sequential (1, 2, ...n) numbering of positions to a "reference" numbering scheme for analysis. This file contains the following columns
 - `sequential_site`: site within glycoprotein ORF, numbered from 1 onwards.
 - `sequential_wt`: wild type amino acid residue at this position
 - `reference_site`: numbering used (conventional in field).
 - `reference_wt`: wild type amino acid residue within reference strain
 - `region`: domain of glycoprotein where position residues, eg. *signal peptide*, *cytoplasmic tail*. 

## Mutation-type classification
[mutation_design_classification.csv](mutation_design_classification.csv) classifies mutations into the different categories of designed mutations. This was generated using `glycoprotein_positions.ipynb` in the scratch_notebook directory.
 - `sequential_site`: site within glycoprotein ORF
 - `amino_acid`: amino acid residue mutation at position
 - `mutation_type`: can be *designed_mutation*, *stop* for introduced stop codons, or *unintended_mutation*.



