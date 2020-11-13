#library(devtools)
#install.packages("RcppArmadillo", repos="https://rcppcore.github.io/drat")
#sudo apt-get install libfftw3-dev libfftw3-doc
#install_github("isglobal-brge/SNPassoc")
library(SNPassoc)

# Next steps may be very time consuming. So they are not executed
#data(HapMap)
#myDat<-setupSNP(HapMap, colSNPs=3:9307, sort = TRUE,
 #info=HapMap.SNPs.pos, sep="")


#resHapMap<-scanWGassociation(group~1, data=myDat, model="log", nperm = 200)

#permTest(resHapMap, 'minimum')


library(tidyverse)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

myY <- read.table('39transformed_WF2.txt', head = TRUE)
myGD <- read_csv('WF2numa0.csv')
myGM <- read_csv('WF2mapa0.csv')
myGM$Chromosome <- sapply(myGM$Chromosome, function(x) as.numeric(gsub("[[:alpha:]]", "", x)))

colnames(myY) <- unname(unlist(sapply(colnames(myY), function(x) str_replace(x, '\\.', ''))))

wGD <- t(myGD[, 2:177])
colnames(wGD) <- unname(unlist(myGD[,1]))

myKI <- read.table("GAPIT.Kin.VanRaden.csv", head = FALSE, sep = ',')

myCV <- read.table("GAPIT.PCA.csv", head = TRUE, sep = ',')


for (file in dir('Original')) {
  df <- read_csv(sprintf('Original/%s', file))
  
  df <- df[order(df$P.value),]
  
  selset <- df$SNP[1:100000]
  
  phen <- unlist(str_split(file, "\\."))[3]
  
  wY <- myY[,c('accession_name', phen)]
    
  oo <- c()
  for (i in 1:500) {
    tGD <- myGD[myGD$taxa %in% selset,]
    wGD2 <- t(tGD[, 2:177])
    colnames(wGD2) <- unname(unlist(tGD[,1]))
    
    wGD2 <- sapply(as.data.frame(wGD2), sample)
    wGD2 <- add_column(as.data.frame(wGD2), rownames(wGD), .before = 1)
    colnames(wGD2)[1] <- 'taxa'
    
    wGM <- myGM[myGM$Name %in% selset,]
    
    myGAPIT <- GAPIT(
      Y=wY,
      GD=wGD2,
      GM=as.matrix(wGM),
      KI=myKI,
      CV=myCV, model = c("Blink"))
    
    df2 <- read_csv(sprintf('GAPIT.Blink.%s.GWAS.Results.csv', phen))
    file.remove(sprintf('GAPIT.Blink.%s.GWAS.Results.csv', phen))
    oo <- c(oo, unlist(unname(df2$P.value[order(df2$P.value)]))[1:100] )
  }
    
   write.table(oo, sprintf('Pout/%s.csv', phen), col.names = FALSE, row.names = FALSE) 
    
  }