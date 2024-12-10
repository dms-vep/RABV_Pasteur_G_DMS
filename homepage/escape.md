---
aside: false
---

# Effects of mutations to rabies G on escape from monoclonal antibodies

This page provides information on the effects of G mutations on escape from a panel of monoclonal antibodies in several ways:

[[toc]]

## Interactive heatmaps

### Heatmaps for individual antibodies

Here are links to interactive heatmaps of the effects of G mutations on escape from neutralization by each monoclonal antibody:

 - [antibody 17C7](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/17C7_mut_effect.html)
 - [antibody CR4098](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/CR4098_mut_effect.html)
 - [antibody RVA122](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/RVA122_mut_effect.html)
 - [antibody RVC58](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/RVC58_mut_effect.html)
 - [antibody CR57](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/CR57_mut_effect.html)
 - [antibody RVC20](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/RVC20_mut_effect.html)
 - [antibody CTB012](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/CTB012_mut_effect.html)
 - [antibody RVC68](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/RVC68_mut_effect.html)

As an example, below is the interactive heatmap showing escape from 17C7; clicking on the above links will open a similar standalone interactive heatmap for each antibody.

<Figure caption="Effects of mutations on escape from neutralization by antibody 17C7">
    <Altair :showShadow="true" :spec-url="'htmls/17C7_mut_effect.html'"></Altair>
</Figure>

The lineplot shows the total escape caused by all amino-acid mutations at each site in G (you can use the interactive options below the plot to show a different site statistics such as mean or max effect of any mutation at a site rather than the sum of mutation effects).
The heatmap shows the escape caused by each mutation, with blue indicating mutations that cause escape (reduce neutralization).
The `x` indicates the parental amino-acid in the Pasteur strain, and dark gray boxes indicate mutations that are too deleterious for cell entry to measure their effect on neutralization, while light gray boxes (which are present for only a few mutations) indicate no reliable measurement was made for that mutation.
You can mouse over mutations in the heatmap for details including the measurement in each of the different experimental replicates.
The mouseovers also report the effect of each mutation on cell entry (`functional effect`).
The `times_seen` parameter in these mousevers represents how many different variants in the pseudovirus library had the mutation, and is an indication of the reliability of the measurement.
Other options at the bottom of the plot allow you to specify additional filters or display options, including whether or not to *floor escape at zero* or show negative escape (mutations that increase neutralization).
You can also filter to only show mutations with cell-entry scores better than some *minimum functional effect*; increasing this slider allows you to exclude mutations that are deleterious for cell entry.
The zoom bar the top can be used to zoom into different regions of protein.

### Lineplots for all antibodies
The above plots show escape for each individual antibody.
[This link](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/all_antibodies_and_cell_entry_faceted.html) provides a set of lineplots that show escape from all of the different antibodies in one interactive plot.

## Escape on structure of G
You can interactively visualize the escape mutations projected onto G's structure using [dms-viz](https://dms-viz.github.io/v0/) at the following links:

The structures at the below links are colored by the total positive escape caused by mutations at each site; you can use the interactive options to plot other site statistics, change the protein representation, etc.
 
### [Prefusion trimer](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fprefusion_ab_escape.json)
Click [here](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fprefusion_ab_escape.json) to visualize escape from all antibodies on the structure of prefusion G (PDB 7u9g).
Initially the lineplot shows escape for all the antibodies overlaid; double-click on an antibody at the toolbar on the left to just show that antibody. Clicking on sites in the lineplot highlights them on the structure.

### [17C7 in complex with G](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fprefusion_17C7_escape.json&bc=%23cfafc7&fi=%257B%2522cell_entry%2522%253A-5%257D)
Click [here](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fprefusion_17C7_escape.json&bc=%23cfafc7&fi=%257B%2522cell_entry%2522%253A-5%257D) to visualize escape from antibody 17C7 mapped onto a structure of the antibody complex with prefusio G (PDB 8a1e).
G can be colored by escape caused by mutations at each site.
Click on sites in the lineplot to show them on the structure of G; the antibody is in pink.

### [RVA22 in complex with G](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fprefusion_RVA122_escape.json&fi=%257B%2522cell_entry%2522%253A-5%257D&bc=%23e1c3e5)
Click [here](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fprefusion_RVA122_escape.json&fi=%257B%2522cell_entry%2522%253A-5%257D&bc=%23e1c3e5) to visualize escape from antibody RVA122 mapped onto a structure of the antibody complex with prefusion G (PDB u9g).
G can be colored by escape caused by mutations at each site.
Click on sites in the lineplot to show them on the structure of G; the antibody is in pink.

### [RVC20 in complex with G pH domain](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fph_domain_RVC20_escape.json&bc=%23c6a9c5&fi=%257B%2522cell_entry%2522%253A-5%257D)
Click [here](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fph_domain_RVC20_escape.json&bc=%23c6a9c5&fi=%257B%2522cell_entry%2522%253A-5%257D) to visualize escape from RVC20 mapped onto a structure of the antibody complex with the pH domain of G (PDB 6tou).
G can be colored by escape caused by mutations at each site.
Click on sites in the lineplot to show them on the structure of G; the antibody is in pink.

### [CR57 in complex with G domain III](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2FdomainIII_CR57_escape.json&bc=%23c9acc5&fi=%257B%2522cell_entry%2522%253A-5%257D)
Click [here](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2FdomainIII_CR57_escape.json&bc=%23c9acc5&fi=%257B%2522cell_entry%2522%253A-5%257D) to visualize escape from CR57 mapped onto a structure of the antibody complex with domain III of G (PDB 8r40).
G can be colored by escape caused by mutations at each site.
Click on sites in the lineplot to show them on the structure of G; the antibody is in pink.

## Numerical values
For the effects of mutations on antibody escape after filtering for only high-confidence measurements, see [this CSV file](https://raw.githubusercontent.com/dms-vep/RABV_Pasteur_G_DMS/refs/heads/main/results/summaries/all_antibodies_and_cell_entry_per_antibody_escape.csv).
If you do not understand the details of the deep mutational scanning experiments and data analysis well enough to understand the data filtering parameters, then we recommend you use that pre-filtered file.

For more detailed values without pre-filtering for high confidence measurements, see the CSV files in [this subdirectory](https://github.com/dms-vep/RABV_Pasteur_G_DMS/tree/main/results/antibody_escape/averages).
Note that someof the mutations in that file may have low `times_seen` (number of variants containing the mutation) or `frac_models` (number of experiments in which the mutation was measured) and so be of lower confidence.
