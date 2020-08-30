library(tidyverse)
tlr <- read_csv('../Data/GWAS data/183x26 wseed.csv')

for (file in dir('../Potential Salt Genes/snpia/')) {
   
snp <- read_csv(sprintf('../Potential Salt Genes/snpia/%s', file), col_types = cols(.default = col_character()))

for (plant in tlr$accession_name) {
  c = 1 
  for (p2 in snp$`JAPONICA NIPPONBARE POSITIONS`)  {
    if (grepl(plant, p2, fixed = TRUE)) {snp[c, 'Id'] = plant}
    c = c + 1 }
  
  c = 1
  for (p2 in snp$ACCESSION) {
    if (grepl(plant, p2, fixed = TRUE)) {snp[c, 'Id'] = plant}
    c <- c + 1 } }

filtered <- snp[-which(is.na(snp$Id)),]
ref <- snp[which(snp$`JAPONICA NIPPONBARE POSITIONS` == "JAPONICA NIPPONBARE ALLELES"),]
ref$Id <- '0'
filtered <- filtered[ - c(which( filtered$`JAPONICA NIPPONBARE POSITIONS` == 'JAO LEUANG::IRGC 65866-1::[IRIS 313-8679]'), which( filtered$`JAPONICA NIPPONBARE POSITIONS` == 'DAMM NOEB THNGORN\'::IRGC 87428-1::[IRIS 313-12152]'), which( filtered$`JAPONICA NIPPONBARE POSITIONS` == "PA PHA EA::IRGC 87751-1::[IRIS 313-12161]")),]

tlr_out = merge(filtered, tlr %>% rename(Id = accession_name), by = 'Id')
tlr_out <- dplyr::bind_rows(tlr_out, ref)

write_csv(tlr_out, sprintf('../Potential Salt Genes/snpiaout/%s', file), na = '') }


