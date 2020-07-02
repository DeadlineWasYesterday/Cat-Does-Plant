#Fisher-Pitman based SNP signifier 

library(tidyverse)
library(coin)

#load data
natol <- read_csv('../Data/GWAS data/natol.csv')

result <- data.frame()
#for everyfile
for (file in dir('../Potential Salt Genes/tsnpout/')) {
  resgen <- data.frame()
  
  snp <- read_csv(sprintf('../Potential Salt Genes/tsnpout/%s', file), col_types = cols(.default = col_character()))
  snp <- snp[-which(is.na(snp$srnm)),]  ##CHANGE HERE
  
  ref <- read_csv(sprintf('../Potential Salt Genes/tsnp/%s', file), col_types = cols(.default = col_character()))
  ref <- ref[1,]
  
  #for every loc
  for (bloc in 7:(length(colnames(snp)) - 8) ) {
    rescol <- data.frame()
    
    #ref base
    orig <- unname(unlist(ref[bloc - 1]))
    
    snips <- na.omit(unique(snp[bloc]))
    for (snip in unname(unlist(snips))) {
      if (nchar(snip) == 1) { std = snip 
      break }}
    snips <- snips[-which(snips == std),]
    
    if (nrow(snips) == 0) {next}
    
    #for every snp
    for (snip in unname(unlist(snips))) {
      if (nchar(snip) == 1) {
      stdlist <- snp$srnm[which(snp[bloc] == std)] ## CHANGE HERE
      sniplist <- snp$srnm[which(snp[bloc] == snip)] ##CHANGE HERE
      
      if(length(stdlist) < 2 || length(sniplist) < 2) {next}
      
      pug <- data.frame(
        val = as.numeric(c(stdlist, sniplist)),
        base = factor(rep(c(std, snip), c(length(stdlist), length(sniplist))))
      )
      
      perm <- pvalue(oneway_test(val ~ base, data = pug, distribution = approximate(nresample = 10000)))
      t <- t.test(as.numeric(stdlist),as.numeric(sniplist))
      
      resnuc <- c(file, paste0(orig,',',std,colnames(snp)[bloc], snip),length(stdlist),length(sniplist),perm,t$p.value)
      resnuc <- unlist(unname(resnuc))
      names(resnuc) <- c('file', 'location', 'lcount', 'rcount', 'perm', 't')
      rescol <- dplyr::bind_rows(rescol, resnuc)
      
      }}
      
      #save filename, blocation and two p values
      resgen <- dplyr::bind_rows(resgen, rescol)
    }
    
    #exit loop and save result
    result <- dplyr::bind_rows(result, resgen)
    
  }
  
#exit loop and save result file

write_csv(result, '../Potential Salt Genes/882x srnm pt.csv')
      
    