#Fisher-Pitman based SNP signifier for all traits

library(tidyverse)
library(coin)

#function that returns perm and t
pnt <- function(std, snip, stdlist, sniplist) {
  pug <- data.frame(val = as.numeric(c(stdlist, sniplist)),
  base = factor(rep(c(std, snip), c(length(stdlist), length(sniplist)))))
  perm <- pvalue(oneway_test(val ~ base, data = pug, distribution = approximate(nresample = 10000)))
  t <- t.test(as.numeric(stdlist),as.numeric(sniplist))
  c(perm, t$p.value) }


for (trt in 1:21) { 

result <- data.frame()
#for everyfile
for (file in dir('../Potential Salt Genes/snpiaout/')) {
  resgen <- data.frame()
  snp <- read_csv(sprintf('../Potential Salt Genes/snpiaout/%s', file), col_types = cols(.default = col_character()))
  
  #for every loc
  for (bloc in 7:(length(colnames(snp)) - 25) ) {
    rescol <- data.frame()
    
    #ref base
    orig <- as.character(snp[which(snp$Id == 0), bloc])
    
    trait <- length(snp) - 21 + trt
    
    sntr <- na.omit(snp[c(bloc, trait)])
    sntr <- sntr[!grepl('/', sntr[1]),] #remove ambigs
    
    snips <- names(table(sntr[1]))[which(table(sntr[1]) > 1)]
    
    if (length(snips) < 2) {next}
    
    refb <- snips[1]
    snips <- snips[-which(snips == refb)]
    
    #for every snp
    for (snip in snips) {
      stdlist <- unlist(unname(sntr[2][which(sntr[1] == refb),])) 
      sniplist <- unlist(unname(sntr[2][which(sntr[1] == snip),]))
        
      pt <- pnt(refb,snip,stdlist,sniplist)
      
      resnuc <- c(file, paste0(orig,',',std,colnames(snp)[bloc], snip),length(stdlist),length(sniplist),pt[1], pt[2])
      resnuc <- unlist(unname(resnuc))
      names(resnuc) <- c('file', 'snp', paste(names(snp[trait]),'lcount'), paste(names(snp[trait]),'rcount'), paste(names(snp[trait]),'p'), paste(names(snp[trait]),'t'))
      rescol <- dplyr::bind_rows(rescol, resnuc)}
  
    resgen <- dplyr::bind_rows(resgen, rescol)
    }
  #exit loop and save result
  result <- dplyr::bind_rows(result, resgen)
}

write_csv(result, sprintf('../Potential Salt Genes/2183x %d pt.csv', trait), na = '')
}

