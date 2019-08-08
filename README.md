# Long Term analysis of influenza virus evolution.

This is the [Nextstrain](https://nextstrain.org) build for seasonal influenza viruses adapted for the publication Le, Neher, Shraiman on ["Phylodynamics of rapidly adapting pathogens: extinction and speciation of a Red Queen."](https://www.biorxiv.org/content/early/2018/10/29/455444).

The results are available as [community builds on nextstrain](https://nextstrain.org/community/neherlab/allflu/h3n2_ha).
The analysis is derived from the pipeline at [nextstrain/seasonal-flu](https://github.com/nextstrain/seasonal-flu). In contrast to the regular nextflu analysis, this analysis operates only on freely available data.

In addition, this repository contains scripts to plot the evolution of the TMRCA for the different lineages.

 * To generate Figure 1, run the Snakemake file. The script `scripts/plot_all_Tmrca.py` generates the saw-tooth graphs in Fig.~1.
 * `script/split_trees.py`: script that splits a global influenza B tree into the sublineages Victoria and Yamagata
 * `script/plot_tmrca.py`: script that uses the output of the `augur` pipeline to produce the tabular files `{lineage}_tmrca_trajectory.dat`. The latter are tab separated files with year in the first and Tmrca in the second columns.

