---
aside: false
---

# Effects of mutations to rabies G on cell entry

This page provides information on the effects of mutations on G-mediated cell entry into 293T cells in several ways:

[[toc]]

## Interactive heatmap

<Figure caption="Effects of mutations on G-mediated cell entry into 293T cells">
    <Altair :showShadow="true" :spec-url="'htmls/HEK293T_entry_func_effects.html'"></Altair>
</Figure>

The heatmap shows the effects of each mutation on cell entry, with negative values (orange) indicating reduced cell entry and positive values (blue) indicating improved cell entry.
The `x` indicates the parental amino-acid in the Pasteur strain, and gray boxes indicate mutations lacking a reliable measurement.
You can mouse over mutations in the heatmap for details including the measurement in each of the different experimental replicates.
The `times_seen` parameter in these mouseovers represents how many different variants in the pseudovirus library had the mutation, and is an indication of the reliability of the measurement.
The lineplot shows a single statistic summarizing the effects of mutations at each site; you can use the options at the bottom of the plot to make this summary statistic the mean, sum, etc of the mutation effects at each site.
Other options at the bottom of the plot allow you to specify additional filters or display options.
The zoom bar the top can be used to zoom into different regions of protein.

Click [here](https://dms-vep.org/RABV_Pasteur_G_DMS/htmls/HEK293T_entry_func_effects.html) for a standalone link to this heatmap.

## Mutation effects on structure of G
You can interactively visualize the effects of mutations on cell entry projected onto G's structure using [dms-viz](https://dms-viz.github.io/v0/) at the following links:

 - [pre-fusion trimer structure PDB 7u9g](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fprefusion_cell_entry.json&sa=true)
 - [low pH extended intermediate monomer PDB 6lgw](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FRABV_Pasteur_G_DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fextended_intermediate_cell_entry.json&sa=true)

The structures at the above links are colored by the mean effect of mutations at each site; you can use the interactive options to plot other site statistics as well.
You can also use the interactive options to select specific sites to show, change the protein representations, etc.

## Numerical values
For the effects of mutations on cell entry after filtering for only high-confidence measurements, see [this CSV file](https://github.com/dms-vep/RABV_Pasteur_G_DMS/blob/main/results/summaries/all_antibodies_and_cell_entry.csv).
If you do not understand the details of the deep mutational scanning experiments and data analysis well enough to understand the data filtering parameters, then we recommend you use that pre-filtered file.

For more detailed values without pre-filtering for high confidence measurements, see [this CSV file](https://github.com/dms-vep/RABV_Pasteur_G_DMS/blob/main/results/func_effects/averages/HEK293T_entry_func_effects.csv); note that a few of the mutations in that file may have low `times_seen` (number of variants containing the mutation) or `n_selections` (number of experiments in which the mutation was measured) and so be of lower confidence.

For the measurements in each of the individual experioments, see the [CSV files in this subdirectory](https://github.com/dms-vep/RABV_Pasteur_G_DMS/tree/main/results/func_effects/by_selection).
