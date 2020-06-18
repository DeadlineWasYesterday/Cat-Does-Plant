

#make result table

result <- data.frame()

#for everyfile
for (file in dir('../Potential Salt Genes/snpiaout/')) {
  
  snp <- read_csv(sprintf('../Potential Salt Genes/snpiaout/%s', file), col_types = cols(.default = col_character()))
  snp$SES <- as.numeric(snp$SES)
  
  ref <- read_csv(sprintf('../Potential Salt Genes/snpia/%s', file), col_types = cols(.default = col_character()))
  ref <- ref[1,]
  
  #split into tolerant and sensitive
  sensitive <- snp[which(snp$SES > 6),]
  tolerant <- snp[which(snp$SES < 6),]
  
  
  #for every column
  resgen <- data.frame()
  for (bloc in 7:(length(colnames(snp)) - 4) ) {
    #matches
    
    senbloc <- sensitive[bloc][!is.na(sensitive[bloc])]
    tolbloc <- tolerant[bloc][!is.na(tolerant[bloc])]
    
    senmat <- sum(senbloc == unlist(unname(ref[bloc - 1])))
    tolmat <- sum(tolbloc == unlist(unname(ref[bloc - 1])))
    senmm <- 0
    tolmm <- 0
    
    
    #mismatches
    for (ltr in senbloc[!(senbloc == unlist(unname(ref[bloc - 1])))]) {
      if (grepl(ltr, unlist(unname(ref[bloc - 1])), fixed = TRUE))
      { senmm <- senmm + 0.25 }
      else { senmm <- senmm + 1}
    }
    
    fqmat <- matrix(c(tolmat, tolmm, senmat, senmm), nrow = 2)
    
    fexact <- fisher.test(fqmat)
    #save filename, blocation and two p values
    rescol <- c(file, colnames(snp)[bloc], fexact['p.value'])
    rescol <- unlist(unname(rescol))
    names(rescol) <- c('file', 'location', 'p')
    resgen <- dplyr::bind_rows(resgen, rescol)
  }
  
  #exit loop and save result
  result <- dplyr::bind_rows(result, resgen)
  
}

#exit loop and save result file
result$p <- as.numeric(result$p)
write.csv(result, '../Potential Salt Genes/882 mmonly split6.csv')



