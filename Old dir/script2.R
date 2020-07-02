#The first three lines import the libraries.

library(tidyverse)
library(seqinr)
#library('ORFik')

#Then the tables are loaded in. 
#Tlr is the table that has tolerance listed for 138 plants. 
#snp is the table you download from the website.

tlr <- read_csv('../Data/Tlr.csv')
snp <- read_csv('../Data/snpi.csv')


#The nested loop is for filtering plants for which tolerance 
#scores are available, from the big snp table. The first for 
#loop loops over the plant ID column in tlr. 
#The second for loop loops over the accession column in the snp table 
#and the third one loops over the first column in the snp table.

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

#The filtered snp table with plants for which tolerance scores 
#are available is saved as filtered. To merge the filtered table 
#with tlr, a third table with some columns removed and some 
#columns renamed is created called tlr_to_merge.

filtered <- filter(snp, snp$Id != 'no')

tlr_to_merge = select(tlr, -X6, -X7, -X8)

tlr_to_merge <- tlr_to_merge %>% rename('Id' = X2)

tlr_out = merge(filtered, tlr_to_merge, by = 'Id')


#The merged table is saved as tlr_out.
#write_csv(tlr_out, '../Data/HKT1;4 SNPi scored.csv')

#tlr_out is loaded in
tlr_out <- read_csv('../Data/HKT1;4 SNPi scored.csv', col_types = cols(.default = col_character()))

#recycled code
#tlr_out$`30735606` <- 'T'
#tlr_out$`30737774` <- 'T'
#tlr_out$`30738142` <- 'T'


#code for using ORFik to find ORFs in the gene sequence.
hf = read.fasta("../Data/OsHKT1;4 genomic.fasta")

str = c2s(hf[[1]])

hf2 = comp(hf[[1]])

str2 = c2s(hf2)

thing <- findORFs(str, startCodon = 'ATG', minimumLength = 400)

thing2 <- findORFs(str, startCodon = 'ATG', minimumLength = 400)


#code for properly specifying positions with uncertain bases
#and missing locations.

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


#recycled code for dismissing rows with missing/uncertain values

#deaths = sum(replace_na(data = tlr_out$kill, replace = FALSE))
#survived <- filter(snp, tlr_out$kill != TRUE)

#write_csv(tlr_out, '../Data/snps.csv')


#merging subsequent rows with insertions
prev_col <- 0
for (col in colnames(tlr_out))
{
  if (grepl('.', col, fixed = TRUE)) {
    tlr_out[prev_col][[1]] <- paste0(tlr_out[prev_col][[1]], tlr_out[col][[1]])
    }
  else { prev_col <- col }
}


#changing chromosome positions to gene positions
c = 1
for (col in colnames(tlr_out))
{
  if (11 <= c  && c <= 226) 
  {
    colnames(tlr_out)[c] <- strtoi(col) - 30734183 #zero indexing
    #print(strtoi(col)-30734183)
  }
  c = c + 1
}

#remove columns previously used for insertions
wtlr <- tlr_out[colnames(tlr_out)[!is.na(colnames(tlr_out))]]

#write snp table with merged insertions and gene positions  
#write_csv(tlr_out, '../Data/HKT1;4 SNPi scored subs.csv')


#load gene sequence
ws = read.fasta("../Data/OsHKT1;4 genRC.fasta")

#execute position specific substitutions with gaps
for (i in 1:(length(rownames(tlr_out)) + 0) )
  {
    tempseq <- read.fasta("../Data/OsHKT1;4 genRC.fasta")
    for (i2 in 11:226)
    {
      tempseq[[1]][strtoi(colnames(tlr_out)[i2]) + 1] <- tlr_out[i,i2][[1]]
  }
    
    ws <- c(ws, tempseq)
}

#ordered naming scheme for fasta files
a = tlr_out$`JAPONICA NIPPONBARE POSITIONS`
a <- c('original', a)

#write a temporary fasta file and load
write.fasta(ws, a, ("../Data/test.fasta"), 
              open = "w", nbchar = 60000, as.string = FALSE)
  
tempseq <- read.fasta("../Data/test.fasta")


#remove gaps and write final file
for (i in 1:length(tempseq)){
tempseq[[i]] <- tempseq[[i]][tempseq[[i]] != "-"] }

write.fasta(tempseq, a, ("../Data/uncurated.fasta"), 
            open = "w", nbchar = 60000, as.string = FALSE)


ws <- read.fasta("../Data/uncurated.fasta")

gs <- lapply(lapply(ws, comp), rev)

ex1 <- lapply(gs, function(x2) x2[85:1168])
ex2 <- lapply(gs, function(x2) x2[3866:4085])
ex3 <- lapply(gs, function(x2) x2[4272:4470])

excds <- lapply(c(1:length(ex1)), function(x) (c(ex1[[x]], ex2[[x]], ex3[[x]])))
       
wt <- lapply(excds, translate)



#wt <- translate(ws[[1]], sens = "R", frame = 0, numcode = 1, ambiguous = TRUE)
#str_length((strsplit(c2s(wt), '*', fixed = TRUE))[[1]])

write.fasta(gs, a, ("../Data/138genes.fasta"), 
            open = "w", nbchar = 60000, as.string = FALSE)
write.fasta(wt, a, ("../Data/138proteins.fasta"), 
            open = "w", nbchar = 60000, as.string = FALSE)


new_tlr <- read_csv('../Data/HKT1;4 SNPi scored.csv', col_types = cols(.default = col_character()))

kill <- c()
for (i in 1:length(rownames(new_tlr))) {
  for (i2 in 11:226) {
    if (!is.na(new_tlr[i,i2])) {
    if (grepl('/', new_tlr[i,i2]) && grepl('-', new_tlr[i,i2])) {
      kill <- c(kill, i)
    }}
  }
}

deaths <- unique(kill)

surv_genes <- gs[-deaths]
surv_proteins <- wt[-deaths]

surv_headers <- a[-deaths]
t129 <- new_tlr[-deaths,]
length(new_tlr[-deaths,])
write_csv(t129, '../Data/129t.csv')

write.fasta(surv_genes, surv_headers, ("../Data/129genes.fasta"), 
            open = "w", nbchar = 60000, as.string = FALSE)
write.fasta(surv_proteins, surv_headers, ("../Data/129proteins.fasta"), 
            open = "w", nbchar = 60000, as.string = FALSE)


