# Deep mutational scanning of the Rabies glycoprotein (G) Pasteur Strain using a barcoded pseudotyped lentiviral platform
Study by Arjun Aditham and Jesse Bloom.


This repo contains data and analyses from deep mutational scanning experiments on the Rabies glycoprotein (G). All experiments were performed on the Pasteur strain of rabies [NC_001542.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001542.1). 


The deep mutational scan only consists of the ectodomain of Rabies G and few sites flanking the ectodomain (positions 18-450). Site 403 is poorly represented. Stop codons were incorporated at alternating positions for sites 18-56 (~20 stop codons). 


For analysis and documentation, navigate to [https://dms-vep.org/RABV_Pasteur_G_DMS/](https://dms-vep.org/RABV_Pasteur_G_DMS/). 

### Data
Input data utilized by the pipeline are located in [./data/](data). 

### Non-pipeline analyses
All other non-pipeline analyses are contained in [./scratch_notebook/](scratch_notebook). The notebooks in this directory are not part of the main pipeline but have been used to generate files used as input for the pipeline.

