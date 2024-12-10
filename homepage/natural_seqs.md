---
aside: false
---

# Antibody escape among natural rabies G sequences

This page provides several ways to identify rabies G sequences that have mutations that are measured in the deep mutational scanning to cause antibody escape:

[[toc]]

## Nextstrain tree
The interactive `nextstrain` tree at [https://nextstrain.org/groups/jbloomlab/dms/rabies-G](https://nextstrain.org/groups/jbloomlab/dms/rabies-G) allows you to examine a phylogeny of rabies G that is colored by the predicted escape of each strain (the sum of the experimentally measured effects of the constituent mutations in the deep mutational scanning).
Note that due to small errors in the values, strains more highly diverged from the Pasteur strain used in the deep mutational scanning will tend to have greater predicted escape.

Use the dropdown at left to choose which antibody to color the tree by.
You can also use the *Genotype* option to color the tree by the amino-acid identity in the *G ectodomain* (those site numbers correspond to the ones used in the deep mutational scanning).

Note that this tree includes lab-adapted as well as natural strains.

## Escape-mutation frequencies
Below is an interactive plot showing the escape caused by each mutation (as measured in the deep mutational scanning) versus its frequency among natural sequences.
You can mouseover points for details about the mutation:

<Figure caption="Escape caused by each mutation versus its frequency among natural sequences">
    <Altair :showShadow="true" :spec-url="'htmls/antibody_escape_vs_natural_diversity.html'"></Altair>
</Figure>
