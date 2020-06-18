library(tidyverse)

tlr <- read_csv('../Data/Tlr.csv', col_types = cols(.default = col_character()))

for (file in dir('../Potential Salt Genes/snpia/')) {
   
snp <- read_csv(sprintf('../Potential Salt Genes/snpia/%s', file), col_types = cols(.default = col_character()))

for (plant in tlr$Plant_Name) {
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
tlr_to_merge = select(tlr, -'992', -'1183', -X8)
#tlr_to_merge <- tlr_to_merge %>% rename('Id' = 'Plant_Name')
colnames(tlr_to_merge)[2] <- 'Id'

tlr_out = merge(filtered, tlr_to_merge, by = 'Id')
write_csv(tlr_out, sprintf('../Potential Salt Genes/snpiaout/%s', file)) }


