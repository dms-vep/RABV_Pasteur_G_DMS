---
layout: home

hero:
  name: "Deep mutational scanning of rabies G (Pasteur strain)"
  tagline: "Pseudovirus deep mutational scanning to measure how rabies G mutations affect cell entry and escape from a panel of monoclonal antibodies"
  image: G_image.png
features:
  - title: Cell entry
    details: Effects of mutations on cell entry
    link: /cell_entry
  - title: Antibody escape
    details: Effects of mutations on escape from monoclonal antibodies
    link: /escape
  - title: Escape in natural sequences
    details: Antibody escape among publicly available rabies G sequences
    link: /natural_seqs
---

## Summary
This website provides easy access to interactive visualizations of the results of [deep mutational scanning of rabies G](https://www.biorxiv.org/content/10.1101/2024.12.17.628970v1).
The links in the boxes above take you to pages that access the data in different forms.

The study that generated these data was performed in the [Bloom lab](https://jbloomlab.org); please see [Aditham et al (2024)](https://www.biorxiv.org/content/10.1101/2024.12.17.628970v1) for full details.

In all the visualizations and data, protein sites are numbered in the standard scheme starting with 1 assigned to the first residue of G's ectodomain (rather than the N-terminal methionine).

The experiments were performed using [pseudovirus deep mutational scanning](https://doi.org/10.1016/j.cell.2023.02.001), which enables the safe study of mutants of viral entry proteins using single-cycle replicative pseudoviruses that are **not** human pathogens.

The GitHub repository with the full analysis and processed data is at [https://github.com/dms-vep/RABV_Pasteur_G_DMS](https://github.com/dms-vep/RABV_Pasteur_G_DMS).
This analysis used [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3/), and the full output of the pipeline (which represents less user-friendly but more extensive documentation) is in the [Appendix](appendix.html){target="_self"}.
