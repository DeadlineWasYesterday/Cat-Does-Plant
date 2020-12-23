#GWAS run on GAPIT

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd('/home/ric/Cat-chan/GWAS_factory/')

#whole population
myY <- read.table('../../Data/39transformed_WF2.txt', head = T)
myG <- read.table('../../Data/WF2.hmp.txt', head = F)

myGAPIT_MLM <- GAPIT(Y=myY, G=myG, PCA.total=3, model=c("Blink", "CMLM"), Multiple_analysis=TRUE, SNP.MAF = 0.05)

#indica subpopulation
myY <- read.table('../../Data/39transformed_WF2ausR7.txt', head = T)
myG <- read.table('../../Data/WF2ausR7.hmp.txt', head = F)

myGAPIT_MLM <- GAPIT(Y=myY, G=myG, PCA.total=3, model=c("Blink", "CMLM"), Multiple_analysis=TRUE, SNP.MAF = 0.1)

#aus subpopulation
myY <- read.table('../../Data/39transformed_WF2indR7.txt', head = T)
myG <- read.table('../../Data/WF2indR7.hmp.txt', head = F)

myGAPIT_MLM <- GAPIT(Y=myY, G=myG, PCA.total=3, model=c("Blink", "CMLM"), Multiple_analysis=TRUE, SNP.MAF = 0.1)