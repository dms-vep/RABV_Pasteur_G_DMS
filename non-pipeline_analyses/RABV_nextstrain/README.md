# Phylogeny analysis for Rabies virus and other Lyssaviruses using Nextstrain

This repo analyzes the current Lyssavirus sequences that have been sequenced by running a [snakemake](https://snakemake.readthedocs.io/) pipeline that downloads sequences based a list of accessions, processes the sequences, and constructs a phylogenetic Nextstrain tree that can be viewed using [Auspice](https://auspice.us/).

Analysis performed by Caleb Carr.

## Genbank accessions

To download the accessions, go to [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and click *Search by virus*. In the *Search by virus name or taxonomy* box, enter *Lyssavirus, taxid:11286* and hit enter. Then click the  *Download* option, select *Accession List* and *Nucleotide* options and hit *Next*. On the next page, select *Download All Records* and hit *Next*. On the next page, select *Accession with version* and click *Download*. Sequences are downloaded from the list of accessions because more information is extracted from the genbank file during the download process. The current accession list was downloaded on June 9, 2024. 

## Snakemake Pipeline

The pipeline can be run automatically using [snakemake](https://snakemake.readthedocs.io/) to run [Snakefile](Snakefile), which reads its configuration from [config.yaml](Configure/config.yml). The results of the automated steps are placed in [Results](Results/).

To run these steps, first build the [conda](https://docs.conda.io/) environment, which installs the necessary programs. First install conda. Then build the environment.

```
conda env create -f environment.yml
```

Then activate the conda environment with:

```
conda activate nextstrain_analysis
```

and then run the pipeline on a computing cluster with [slurm](https://slurm.schedmd.com/documentation.html), which uses the configuration specified in [cluster.yml](cluster.yml):

```
sbatch run_snakemake_cluster.bash
```

## Organization of this repo

- [`Configure`](Configure/): Contains the files needed to configure the pipeline including the [config.yaml](Configure/config.yml) and [Input_Data](Configure/Input_Data/) which contains list of accessions, reference genomes, and outgroup references. 
- [`Rules`](Rules/): Contains the snakemake rules used to run the pipeline.
- [`Scripts`](Scripts/): Contains the custom [python](https://www.python.org/) scripts used for part of the analysis.
- [`Results`](Results/): Contains the automatically generated results from the pipeline analysis.