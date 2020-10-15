#WF2Miss25MAF10GPBLINK

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd('/home/ric/Cat-chan/GWAS_factory/WF2Miss25MAF10GPBLINK')

myY <- read.table('../../Data/38ProcessedTraits.txt', head = T)
myG <- read.table('../../Data/WF2Miss25.hmp.txt', head = F)

myGAPIT_MLM <- GAPIT( Y=myY, G=myG, PCA.total=3, model=c("Blink"), Multiple_analysis=TRUE, Model.selection = TRUE, SNP.MAF = 0.1, SNP.FDR = 0.01)