library(tidyverse)

#Make sure all tables have snps
for (file in list.files('../Potential Salt Genes/tsnp/')) {
  snp <- read_csv(sprintf('../Potential Salt Genes/tsnp/%s', file), col_types = cols(.default = col_character()))
    if (! grepl('.', gsub(", ","",toString(colnames(snp))), fixed = TRUE)) {
      print(file)
    }
} #Only one gene reported to have no indels CHR3-11774086-11778008
#and CHR5-7474468-7477040 has no insertions


#find files
rapt <- read_csv('../Potential Salt Genes/Base/24200×27 first hit genelist.csv')

for (file in list.files('../Potential Salt Genes/snpia/')) {
  chr <- unlist(str_split(file, '-'))[2]
  wps <- as.integer(unlist(str_split(file, '-'))[3])
  end <- as.integer(unlist(str_split(file, '-'))[4])
  
  l1 <- rapt[which(rapt$Chromosome == chr),]
  l2 <- l1[which(l1$`With Promoter Start` >= wps),]
  l3 <- l2[which(l2$End <= end),]
  check <- l3$Gene[1]
  
  rapt$Status[which(rapt$Gene == check)] <- 'Done.'
   }


#rename
rapt <- rapt[rapt$Status == 'Done.',]
rapt$wGName <- unlist(lapply(rapt$`Oryzabase Gene Symbol Synonym(s)`, str_extract, pattern = '([A-Za-z0-9])+'))
rapt$Identifier <- paste0(rapt$Chromosome, '-', rapt$`With Promoter Start`, '-', rapt$End)

#remove duplicates
#View(rapt[which(duplicated(rapt$Identifier)),])

rapt <- rapt[-which(duplicated(rapt$Identifier)),] #kills the list when no duplicates
rapt$wGName[which(is.na(rapt$wGName))] <- paste0('NoName', c(1:length(which(is.na(rapt$wGName)))))

rapt$filename <- paste0(rapt$Gene, '-', rapt$Identifier, '-', rapt$wGName)

#check duplicate files and delete
files <- list.files('../Potential Salt Genes/snpia/')
filei <- (unlist(lapply(files, str_extract, pattern = '(-.*-)')))
fdf <- data.frame(files, filei)

for (hit in fdf$files[duplicated(fdf$filei)]) {
  file.remove(sprintf('../Potential Salt Genes/snpia/%s', hit)) }


rapt$MSUfilename <- ''
#second filename column for msu files
#needs to be run after duplicates are removed
for (file in list.files('../Potential Salt Genes/snpia/')) {
  chr <- unlist(str_split(file, '-'))[2]
  wps <- as.integer(unlist(str_split(file, '-'))[3])
  end <- as.integer(unlist(str_split(file, '-'))[4])
  
  l1 <- rapt[which(rapt$Chromosome == chr),]
  l2 <- l1[which(l1$`With Promoter Start` >= wps),]
  l3 <- l2[which(l2$End <= end),]
  check <- l3$Gene[1]
  
  rapt$MSUfilename[which(rapt$Gene == check)] <- file }                                        


#rename
for (file in rapt$MSUfilename) {
  file.rename(from=sprintf('../Potential Salt Genes/snpia/%s', file),to=sprintf('../Potential Salt Genes/snpia/%s.csv', rapt$filename[which(rapt$MSUfilename == file)])) }


#order 

rapt <- rapt[order(rapt$Gene),]

#compile columns into bfile
bfl <- read_csv("../Potential Salt Genes/snpiaout/Os04t0473900-CHR4-23696662-23707348-DWA1.csv", col_types = cols(.default = col_character()))
bfl <- bfl[order(bfl['Id']),]
bfl <- bfl[c(1:5, (length(bfl)-24):(length(bfl)))]
for (file in list.files('../Potential Salt Genes/snpiaout/')) {
  snp <- read_csv(sprintf('../Potential Salt Genes/snpiaout/%s', file), col_types = cols(.default = col_character()))
  snp <- snp[order(snp['Id']),]
  if (strtrim(file, 4) == 'Os02') {
  for (cn in 7:(length(colnames(snp)) - 25)) {
    snips <- unlist(unique(na.omit(snp[cn])))
    if (length(snips[!grepl('/', snips)]) > 1){
      if (sum(colnames(bfl) == names(snp[cn])) != 0) {next}
      bfl <- dplyr::bind_cols(bfl, snp[cn])
    }
    
  }}
  if (strtrim(file, 4) == 'Os03') {break}
}

write_csv(bfl, '../Potential Salt Genes/bfch2.csv', na = '')




#check missing rows
for (file in list.files('../Potential Salt Genes/snpiaout/')) {
  snp <- read_csv(sprintf('../Potential Salt Genes/snpiaout/%s', file), col_types = cols(.default = col_character()))
  if (nrow(snp) != 175) {print(file)}
}
