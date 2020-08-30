#Na tolerance

library(tidyverse)

subgrp <- read_csv('../Data/GWAS data/bothgrps.csv')

csesgrp <- read_csv('../Data/GWAS data/csesgrp.csv')

natol <- merge(select(subgrp, Id = accession_name,srnm,ssnm,rwm,srwm,swm,sswm), 
               select(csesgrp, Id = accession_name,pcaloss,pcbloss,SESm),
               by = 'Id')
natol['prwloss'] <- ( natol$rwm - natol$srwm ) / natol$rwm * 100
natol['pswloss'] <- ( natol$swm - natol$sswm ) / natol$swm * 100

natol <- select(natol,Id,srnm,ssnm,SESm,pcaloss,pcbloss,prwloss,pswloss)

#write_csv(natol, '../Data/GWAS data/natol.csv')
natol <- read_csv('../Data/GWAS data/natol.csv')

##

for (file in dir('../Potential Salt Genes/tsnp/')) {
  snp <- read_csv(sprintf('../Potential Salt Genes/tsnp/%s', file), col_types = cols(.default = col_character()))
  for (plant in natol$Id) {
    c = 1 
    for (p2 in snp$`JAPONICA NIPPONBARE POSITIONS`)  {
      if (grepl(plant, p2, fixed = TRUE)) {snp[c, 'Id'] = plant}
      c = c + 1 }
    c = 1
    for (p2 in snp$ACCESSION) {
      if (grepl(plant, p2, fixed = TRUE)) {snp[c, 'Id'] = plant}
      c <- c + 1 } }
  
  filtered <- snp[-which(is.na(snp$Id)),]
  filtered <- filtered[ - c(which( filtered$`JAPONICA NIPPONBARE POSITIONS` == 'JAO LEUANG::IRGC 65866-1::[IRIS 313-8679]'), which( filtered$`JAPONICA NIPPONBARE POSITIONS` == 'DAMM NOEB THNGORN\'::IRGC 87428-1::[IRIS 313-12152]')),]
  
  tlr_out <- merge(filtered, natol, by = 'Id')
  
  write_csv(tlr_out, sprintf('../Potential Salt Genes/tsnpout/%s', file)) }


####


result <- data.frame()
#for everyfile
for (file in dir('../Potential Salt Genes/tsnpout/')) {
  resgen <- data.frame()
  
  snp <- read_csv(sprintf('../Potential Salt Genes/tsnpout/%s', file), col_types = cols(.default = col_character()))
  #snp$pcaloss <- as.numeric(snp$pcaloss)
  
  ref <- read_csv(sprintf('../Potential Salt Genes/tsnp/%s', file), col_types = cols(.default = col_character()))
  ref <- ref[1,]
  
  #split into tolerant and sensitive
  sensitive <- snp[which(snp$pcbloss > 35),]
  sensitive <- dplyr::bind_rows(sensitive, ref)
  tolerant <- snp[which(snp$pcbloss < 0),]
  
  #for every column
  for (bloc in 7:(length(colnames(snp)) - 8) ) {
    rescol <- data.frame()
    
    #ref base
    orig <- unname(unlist(ref[bloc - 1]))
    
    #count whole letters
    asen <- sum(sensitive[bloc][!is.na(sensitive[bloc])] == 'A')
    atol <- sum(tolerant[bloc][!is.na(tolerant[bloc])] == 'A')
    
    csen <- sum(sensitive[bloc][!is.na(sensitive[bloc])] == 'C')
    ctol <- sum(tolerant[bloc][!is.na(tolerant[bloc])] == 'C')
    
    tsen <- sum(sensitive[bloc][!is.na(sensitive[bloc])] == 'T')
    ttol <- sum(tolerant[bloc][!is.na(tolerant[bloc])] == 'T')
    
    gsen <- sum(sensitive[bloc][!is.na(sensitive[bloc])] == 'G')
    gtol <- sum(tolerant[bloc][!is.na(tolerant[bloc])] == 'G')
    
    gapsen <- sum(sensitive[bloc][!is.na(sensitive[bloc])] == '-')
    gaptol <- sum(tolerant[bloc][!is.na(tolerant[bloc])] == '-')
    
    #count half letters
    wc <- unlist(unname(sensitive[bloc]))
    ixr <- lapply(wc, grepl, pattern = '/', fixed = TRUE)
    ixr <- as.logical(as.numeric(unlist(ixr)))
    hls <- wc[ixr]
    lts <- unlist(strsplit(hls, '/'))
    hcasen <- sum(lts == 'A')
    hccsen <- sum(lts == 'C')
    hctsen <- sum(lts == 'T')
    hcgsen <- sum(lts == 'G')
    hcgapsen <- sum(lts == '-')
    
    wc <- unlist(unname(tolerant[bloc]))
    ixr <- lapply(wc, grepl, pattern = '/', fixed = TRUE)
    ixr <- as.logical(as.numeric(unlist(ixr)))
    hls <- wc[ixr]
    lts <- unlist(strsplit(hls, '/'))
    hcatol <- sum(lts == 'A')
    hcctol <- sum(lts == 'C')
    hcttol <- sum(lts == 'T')
    hcgtol <- sum(lts == 'G')
    hcgaptol <- sum(lts == '-')
    
    #count total letters
    totasen <- asen + hcasen/2
    totcsen <- csen + hccsen/2
    totgsen <- gsen + hcgsen/2
    tottsen <- tsen + hctsen/2
    totgapsen <- gapsen + hcgapsen/2
    totatol <- atol + hcatol/2
    totctol <- ctol + hcctol/2
    totgtol <- gtol + hcgtol/2
    totttol <- ttol + hcttol/2
    totgaptol <- gaptol + hcgaptol
    
    #make freq table
    fqmat <- matrix(c(totasen,totcsen,totgsen,tottsen,totgapsen,totatol,totctol,totgtol,totttol,totgaptol), nrow = 5)
    fqtab <- data.frame(fqmat)
    rownames(fqtab) <- c('A', 'C', 'G', 'T', '-')
    fqtab <- fqtab[rowSums(fqtab) != 0,]
    
    #for every snp
    for (nuc in rownames(fqtab)[rownames(fqtab) != orig]) { 
      resnuc <- data.frame()
      square <- t(matrix(c(as.numeric(fqtab[orig,]), as.numeric(fqtab[nuc,])), nrow = 2))
      
      #exception handling
      if (sum(square) == 0) {square <-(matrix(c(1,1,1,1),nrow=2))}
      
      
      #Calculate Fischer's exact
      fexact <- fisher.test(square)
      
      resnuc <- c(file, paste0(orig, colnames(snp)[bloc], nuc),  fexact['p.value'])
      resnuc <- unlist(unname(resnuc))
      names(resnuc) <- c('file', 'location', 'p')
      rescol <- dplyr::bind_rows(rescol, resnuc)
      
    }
    #save filename, blocation and two p values
    resgen <- dplyr::bind_rows(resgen, rescol)
  }
  
  #exit loop and save result
  result <- dplyr::bind_rows(result, resgen)
  
}

#exit loop and save result file
result$p <- as.numeric(result$p)
write.csv(result, '../Potential Salt Genes/882x pcbl l0m35.csv')

