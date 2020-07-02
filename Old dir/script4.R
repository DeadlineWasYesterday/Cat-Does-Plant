library(ape)
library(phangorn)
library(seqinr)

dat <- read.dna("../Data/129g_saligned.fasta", format="fasta")
dat_phyDat <- phyDat(dat, type = "DNA", levels = NULL)

dna_dist <- dist.ml(dat, model="JC69")

dat_UPGMA <- upgma(dna_dist)
dat_NJ  <- NJ(dna_dist)

plot.phylo(dat_UPGMA, cex = 0.15)

plot.phylo(dat_NJ, cex = 0.15)





