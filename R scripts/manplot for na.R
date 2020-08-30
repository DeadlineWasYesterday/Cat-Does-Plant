library(tidyverse)

dat <- read_csv('../Data/GWAS data/wf.csv')

accgrp <- group_by(dat, accession_name, Treatment)

accgrp <- accgrp[which(accgrp$`Total replicate` >= 3),]

rsn <- summarise(accgrp, rnm = mean(`Root Na+`), rns = sd(`Root Na+`),
        snm = mean(`Shoot Na+`), sns = sd(`Shoot Na+`),
        cam = mean(`Chlorophyll A`), cas = sd(`Chlorophyll B`),
        cbm = mean(`Chlorophyll B`), cbs = sd(`Chlorophyll B`),
        sesm = mean(`SES`), sess = sd(`SES`),
        rwm = mean(`Root weight`), rws = sd(`Root weight`),
        swm = mean(`Shoot weight`), sws = sd(`Shoot weight`),
        rlm = mean(`Root length`), rls = sd(`Root length`),
        slm = mean(`Shoot length`), sls = mean(`Shoot length`))

congrp <- rsn[which(rsn$Treatment == 'Control'),]
#write_csv(congrp, '../Data/GWAS data/congrp.csv')
strgrp <- rsn[which(rsn$Treatment == 'Stress'),]
#write_csv(strgrp, '../Data/GWAS data/strgrp.csv')

subgrp <- merge(select(congrp, -Treatment, -sesm, -sess) ,select(strgrp, accession_name,
srnm = rnm, srns = rns, ssnm = snm, ssns = sns, scam = cam, scas = cas,
scbm = cbm, scbs = cbs, SESm = sesm, SESs = sess, srwm = rwm, srws = rws,
sswm = swm, ssws = sws, srlm = rlm, srls = rls, sslm = slm, ssls = sls,
), by= 'accession_name')

#write.csv(subgrp, '../Data/GWAS data/bothgrps.csv')

subgrp['pcaloss'] <- ( subgrp$cam - subgrp$scam ) / subgrp$cam
subgrp['pcaloss'] <- ( subgrp$cam - subgrp$scam ) / subgrp$cam * 100
subgrp['pcbloss'] <- ( subgrp$cbm - subgrp$scbm ) / subgrp$cbm * 100

subgrp['l2rn'] <- -log2(subgrp$rnm)
subgrp['l2srn'] <- -log2(subgrp$srnm)
subgrp['l2sn'] <- -log2(subgrp$snm)
subgrp['l2ssn'] <- -log2(subgrp$ssnm)

subgrp['root na fold increase'] <- subgrp$l2rn - subgrp$l2srn
subgrp['shoot na fold increase'] <- subgrp$l2sn - subgrp$l2ssn

csesgrp <- select(subgrp, accession_name,cam,cbm,scam,scbm,SESm, pcaloss, pcbloss)
#write_csv(csesgrp, '../Data/GWAS data/csesgrp.csv')

#correlation matrices
congrp <- read_csv('../Data/GWAS data/congrp.csv')
strgrp <- read_csv('../Data/GWAS data/strgrp.csv')
subgrp <- read_csv('../Data/GWAS data/bothgrps.csv')

nonacongrp <- (na.omit(select(congrp, -sesm, -sess)))

cormat <- round(cor(nonasubgrp[1:17 * 2]),2)


nonasubgrp <- na.omit(subgrp)

library(reshape2)
melted_cormat <- melt(cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "grey")+
  scale_fill_gradient2(low = "blue", high = "green", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson's\nCorrelation") +
  theme_minimal()+ 
  coord_fixed()








###########################

#ses vs shoot na
temp <- select(rsn, sesm, sess, 
               snm, sns)
temp['sesv'] <- temp$sess/ temp$sesm
temp['snv'] <- temp$sns / temp$snm
temp <- temp[which(temp$sesv < 0.2),]
temp <- temp[which(temp$snv < 0.3),]
cor.test(temp$sesm, temp$snm)

#ses vs chl a
temp <- select(rsn, sesm, sess, 
               cam, cas)
temp['sesv'] <- temp$sess/ temp$sesm
temp['cav'] <- temp$cas / temp$cam
temp <- temp[which(temp$sesv < 0.35),]
temp <- temp[which(temp$cav < 0.9),]
cor.test(temp$sesm, temp$cam)




rsnc <- rsn[which(rsn$Treatment == 'Control'),]
rsns <- rsn[which(rsn$Treatment == 'Stress'),]

#rsnc <- rsnc[!is.na(rsnc$rnm),]
#rsns <- rsns[!is.na(rsns$rnm),]

crsnc <- rsnc[rsnc$accession_name %in% rsns$accession_name,]

crsns <- rsns[rsns$accession_name %in% rsnc$accession_name,]

sub <- crsns[order(crsns$accession_name),][c(3,5,7,9,13,15,17,19)] - crsnc[order(crsnc$accession_name),][c(3,5,7,9,13,15,17,19)]


sub <- bind_cols(crsns[order(crsns$accession_name),][1], sub)
colnames(sub)[1] <- 'Id'
#end of nasub


rsm <-  summarise(accgrp, snm = mean(`Shoot Na+`), sns = sd(`Shoot Na+`),
                  cam = mean(`Chlorophyll A`), cas = sd(`Chlorophyll B`),
                  cbm = mean(`Chlorophyll B`), cbs = sd(`Chlorophyll B`),
                  sesm = mean(`SES`), sess = sd(`SES`))

for_ses <- rsm[which(rsm$sess < 2),]




tlr <- read_csv('../Data/Tlr.csv', col_types = cols(.default = col_character()))
tlr <- tlr %>% rename('Id' = 'Plant_Name')
tlr <- select(tlr, -'992', -'1183', -X8)

for (file in dir('../Potential Salt Genes/tsnp/')) {
  
  snp <- read_csv(sprintf('../Potential Salt Genes/tsnp/%s', file), col_types = cols(.default = col_character()))
  
  for (plant in sub$Id) {
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

  
  
  tlr_out <- merge(filtered, sub, by = 'Id')
  tlr_out <- merge(tlr_out, tlr, by = 'Id', all = TRUE)
  tlr_out <- tlr_out[-which(is.na(tlr_out$snm)),]
  
  write_csv(tlr_out, sprintf('../Potential Salt Genes/tsnpout/%s', file)) }



##################################################################

result <- data.frame()

#for everyfile
for (file in dir('../Potential Salt Genes/tsnpout/')) {
  resgen <- data.frame()
  
  snp <- read_csv(sprintf('../Potential Salt Genes/tsnpout/%s', file), col_types = cols(.default = col_character()))
  snp$snm <- as.numeric(snp$snm)
  
  ref <- read_csv(sprintf('../Potential Salt Genes/tsnp/%s', file), col_types = cols(.default = col_character()))
  ref <- ref[1,]
  
  #split into tolerant and sensitive
  sensitive <- snp[which(snp$snm > 0.15),]
  sensitive <- dplyr::bind_rows(sensitive, ref)
  tolerant <- snp[which(snp$snm < 0.07),]
  
  #for every column
  for (bloc in 7:(length(colnames(snp)) - 6) ) {
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
write.csv(result, '../Potential Salt Genes/882x snm 07x015.csv')

