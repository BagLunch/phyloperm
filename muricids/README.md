# muricids
Data from muricid case study along with R code for analysing it.
Data from Pappalardo, Paula, Enrique Rodriguez-Serrano, and Miriam Fernandez. "Correlated evolution between mode of larval development and habitat in muricid gastropods." PloS one 9.4 (2014)

- Table_S1_Pappalardo_et_al_2014.csv
Self-explanatory. Used in logistic PGLS, but not the phylogenetic permutation test.

- muricidae.dated.tre
Phylogeny time-scaled with treepl (see "muricidae_treepl")

- muricids.R
Code for running all analyses and generating all figures. The function "phylo.permute" is in a different folder.

- rangesmuri_phylogeneticAnalysis.csv
Geographic range data used in generating phylogenetic permutations and analyzing ecogeographic data "as such."

- muricidae_treepl
Contains everything needed to time-scale the phylogeny from Pappalardo et al.
