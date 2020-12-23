#Heritability estimates from kinship

library(heritability)
library(tidyverse)

myY <- read_csv('39traitsord.csv')
gv <- pull(myY[,1])

k <- read_csv('GAPIT.Kin.VanRaden.csv', col_names = FALSE)

tn <- k[,1]

k <- as.matrix(k[,2:177])
rownames(k) <- pull(tn)
colnames(k) <- pull(tn)

ov <- c()
for (i in 2:40) {
dv <- pull(myY[,i])

out <- marker_h2(data.vector = dv, geno.vector = gv, K = k)

ov <- c(ov, out$h2)
}

odf <- dplyr::bind_cols(colnames(myY)[2:40], ov)
colnames(odf) <- c('Trait', 'Heritability')

write_csv(odf, 'heritability.csv')
