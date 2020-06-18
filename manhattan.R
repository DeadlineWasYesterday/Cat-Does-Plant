

#make result table

result <- data.frame()

#for everyfile
for (file in dir('../Potential Salt Genes/snpiaout/')) {
  
  resgen <- data.frame()
  snp <- read_csv(sprintf('../Potential Salt Genes/snpiaout/%s', file), col_types = cols(.default = col_character()))
  snp$SES <- as.numeric(snp$SES)
  
  #split into tolerant and sensitive
  sensitive <- snp[which(snp$SES > 8),]
  tolerant <- snp[which(snp$SES < 8),]

  #for every column
  for (bloc in 7:(length(colnames(snp)) - 4) ) {
    
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
    
    #supposedly faulty. omitted.
    #fqmat <- matrix(c(totasen,totatol,totcsen,totctol,totgsen,totgtol,tottsen,totttol,totgapsen,totgaptol), nrow = 5)
    #remove zero rows. not mandatory
    #fqmat[rowSums(fqmat) != 0,]
    fqmat <- matrix(c(totasen,totcsen,totgsen,tottsen,totgapsen,totatol,totctol,totgtol,totttol,totgaptol), nrow = 5)
    
    #Calculate Fischer's exact
    
    #exception handling
    if (sum(fqmat) == 0) {fqmat <-(matrix(c(1,1,1,1),nrow=2))}
    
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
write.csv(result, '../Potential Salt Genes/882 at tlr split8.csv')





## Taking only integer SES

#make result table

result <- data.frame()

#for everyfile
for (file in dir('../Potential Salt Genes/snpiaout/')) {
  
  resgen <- data.frame()
  snp <- read_csv(sprintf('../Potential Salt Genes/snpiaout/%s', file), col_types = cols(.default = col_character()))
  snp$SES <- as.numeric(snp$SES)
  
  #remove fracs
  snp <- snp[which(snp$SES - as.integer(snp$SES) == 0),]
  
  #split into tolerant and sensitive
  sensitive <- snp[which(snp$SES > 4),]
  tolerant <- snp[which(snp$SES < 4),]
  
  #for every column
  for (bloc in 7:(length(colnames(snp)) - 4) ) {
    
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
    
    #supposedly faulty. omitted.
    #fqmat <- matrix(c(totasen,totatol,totcsen,totctol,totgsen,totgtol,tottsen,totttol,totgapsen,totgaptol), nrow = 5)
    #remove zero rows. not mandatory
    #fqmat[rowSums(fqmat) != 0,]
    fqmat <- matrix(c(totasen,totcsen,totgsen,tottsen,totgapsen,totatol,totctol,totgtol,totttol,totgaptol), nrow = 5)
    
    #Calculate Fischer's exact
    
    #exception handling
    if (sum(fqmat) == 0) {fqmat <-(matrix(c(1,1,1,1),nrow=2))}
    
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
write.csv(result, '../Potential Salt Genes/882 nofrac attlr split4.csv')
