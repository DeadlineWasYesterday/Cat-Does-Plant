library(tidyverse)

natol <- read_csv('../Data/GWAS data/natol.csv')

shapiro.test(natol$srnm)
shapiro.test(natol$ssnm)
shapiro.test(natol$SESm)
shapiro.test(natol$pcaloss)
shapiro.test(natol$pcbloss)
shapiro.test(natol$prwloss)
shapiro.test(natol$pswloss)

pcal0g25 <- read_csv('../Potential Salt Genes/882/882x pcbl l0m35.csv')
pcal0g25 <- pcal0g25[which(pcal0g25$p < 0.15),]

files <- unique(pcal0g25$file)

result <- data.frame()
for (file in files) {
  wp <- pcal0g25[which(pcal0g25$file == file),]
  snp <- read_csv(sprintf('../Potential Salt Genes/tsnpout/%s', file), col_types = cols(.default = col_character()))
  
  for (loc in wp$location) {
    lref <- str_extract(loc, '^([A-Z-])')
    lnum <- str_extract(loc, '[0-9]+')
    lpol <- str_extract(loc, '([A-Z-])$')
      
    refset <- snp[which(snp[lnum] == lref),]
    polset <- snp[which(snp[lnum] == lpol),]
    if (nrow(refset) < 2 || nrow(polset) < 2){ next }
    wp[which(wp$location == loc), 'lenrefset'] <- length(rownames(refset))
    wp[which(wp$location == loc), 'lenpolset'] <- length(rownames(polset))
    
    #omitna 
    if (length(na.omit(polset$ssnm)) > 1 && length(na.omit(refset$ssnm)) > 1) {
    wp[which(wp$location == loc), 'tofsnm'] <- t.test(as.numeric(refset$ssnm), as.numeric(polset$ssnm))$p.value }
    if (length(na.omit(polset$pcbloss)) > 1 && length(na.omit(refset$pcbloss)) > 1) {
    wp[which(wp$location == loc), 'tofpcal'] <- t.test(as.numeric(refset$pcaloss), as.numeric(polset$pcaloss))$p.value
    wp[which(wp$location == loc), 'tofpcbl'] <- t.test(as.numeric(refset$pcbloss), as.numeric(polset$pcbloss))$p.value }
    
    wp[which(wp$location == loc), 'tofSES'] <- t.test(as.numeric(refset$SESm), as.numeric(polset$SESm))$p.value
    wp[which(wp$location == loc), 'tofprwl'] <- t.test(as.numeric(refset$prwloss), as.numeric(polset$prwloss))$p.value
    wp[which(wp$location == loc), 'tofpswl'] <- t.test(as.numeric(refset$pswloss), as.numeric(polset$pswloss))$p.value
  }
  
  result <- dplyr::bind_rows(result, wp)
  
}

result[7:12][is.na(result[7:12])] <- 1
result$score <- -log10(result$tofsnm) + -log10(result$tofpcal) + -log10(result$tofpcbl) + -log10(result$tofSES) + -log10(result$tofprwl) + -log10(result$tofpswl)

write_csv(result, '../Potential Salt Genes/882x pcbl l0m35 with t.csv')
