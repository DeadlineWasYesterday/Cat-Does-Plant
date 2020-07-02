library(tidyverse)
library(seqinr)
#library('ORFik')

tlr <- read_csv('../Data/Tlr.csv')
snp <- read_csv('../Data/snpi.csv')

for (plant in tlr$X2)
{
  c = 1 
  
  for (p2 in snp$`JAPONICA NIPPONBARE POSITIONS`)
  {
    
    if (grepl(plant, p2, fixed = TRUE))
      
    {snp[c, 'Id'] = plant}
    c = c + 1
    
  }
  
  c = 1
  for (p2 in snp$ACCESSION)
  {
    if (grepl(plant, p2, fixed = TRUE))
      
    {snp[c, 'Id'] = plant}
    
    c <- c + 1
  }
    
}

filtered <- filter(snp, snp$Id != 'no')

tlr_to_merge = select(tlr, -X6, -X7, -X8)

tlr_to_merge <- tlr_to_merge %>% rename('Id' = X2)

tlr_out = merge(filtered, tlr_to_merge, by = 'Id')

#write_csv(tlr_out, '../Data/HKT1;4 SNPi scored.csv')
tlr_out <- read_csv('../Data/HKT1;4 SNPi scored.csv', col_types = cols(.default = col_character()))

#tlr_out$`30735606` <- 'T'
#tlr_out$`30737774` <- 'T'
#tlr_out$`30738142` <- 'T'

hf = read.fasta("../Data/OsHKT1;4 genomic.fasta")

str = c2s(hf[[1]])

hf2 = comp(hf[[1]])

str2 = c2s(hf2)

thing <- findORFs(str, startCodon = 'ATG', minimumLength = 400)

thing2 <- findORFs(str, startCodon = 'ATG', minimumLength = 400)


for (col in 11:226) {
  for (i in 1:(length(rownames(tlr_out)))) {
    if (grepl('/', tlr_out[i, col], fixed = TRUE)){
      if (grepl('-', tlr_out[i, col], fixed = TRUE)) {
        tlr_out[i, col] <- 'N' }
      else if (grepl('A', tlr_out[i, col], fixed = TRUE) && 
          (grepl('T', tlr_out[i, col], fixed = TRUE))) {
        tlr_out[i, col] <- 'W' }
      else if (grepl('C', tlr_out[i, col], fixed = TRUE) && 
               (grepl('G', tlr_out[i, col], fixed = TRUE))) {
        tlr_out[i, col] <- 'S' }
      else if (grepl('A', tlr_out[i, col], fixed = TRUE) && 
               (grepl('C', tlr_out[i, col], fixed = TRUE))) {
        tlr_out[i, col] <- 'K' }
      else if (grepl('G', tlr_out[i, col], fixed = TRUE) && 
               (grepl('T', tlr_out[i, col], fixed = TRUE))) {
        tlr_out[i, col] <- 'M' }
      else if (grepl('A', tlr_out[i, col], fixed = TRUE) && 
               (grepl('G', tlr_out[i, col], fixed = TRUE))) {
        tlr_out[i, col] <- 'Y' }
      else if (grepl('C', tlr_out[i, col], fixed = TRUE) && 
                 (grepl('T', tlr_out[i, col], fixed = TRUE))) {
        tlr_out[i, col] <- 'R' }
    }
    else if (is.na(tlr_out[i, col])) {
      tlr_out[i, col] <- 'N' }
    }}

deaths = sum(replace_na(data = tlr_out$kill, replace = FALSE))


survived <- filter(snp, tlr_out$kill != TRUE)

write_csv(tlr_out, '../Data/snps.csv')

prev_col <- 0
for (col in colnames(tlr_out))
{
  
  if (grepl('.', col, fixed = TRUE)) {
    tlr_out[prev_col][[1]] <- paste0(tlr_out[prev_col][[1]], tlr_out[col][[1]])
    
  }
  else { prev_col <- col }

}



c = 1
for (col in colnames(tlr_out))
{
  if (11 <= c  && c <= 71) 
  {
    colnames(tlr_out)[c] <- strtoi(col) - 30734183 #zero indexing
    #print(strtoi(col)-30734183)
  }
  c = c + 1
}
  
#write_csv(tlr_out, '../Data/HKT1;4 SNP scored subs.csv')

  ws = read.fasta("../Data/OsHKT1;4 genRC.fasta")
  
  for (i in 1:(length(rownames(tlr_out)) + 0) )
  {
    tempseq <- read.fasta("../Data/OsHKT1;4 genRC.fasta")
    for (i2 in 11:226)
    {
      if (is.na(tlr_out[i,i2]) )
      {
        
      }
      else if (tlr_out[i,i2] == 'A' || tlr_out[i,i2] == 'T'
          || tlr_out[i,i2] == 'G' || tlr_out[i,i2] == 'C'
          || tlr_out[i,i2] == '-')
      {
    
    tempseq[[1]][strtoi(colnames(tlr_out)[i2]) + 1] <- tlr_out[i,i2][[1]]
  }
    }
    ws <- c(ws, tempseq)
  }
  
  a = tlr_out$`JAPONICA NIPPONBARE POSITIONS`
  a <- c('original', a)
  
  write.fasta(ws, a, ("../Data/test.fasta"), 
              open = "w", nbchar = 60000, as.string = FALSE)
