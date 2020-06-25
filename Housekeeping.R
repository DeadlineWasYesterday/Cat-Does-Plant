library(tidyverse)

#Make sure all tables have snps
for (file in list.files('../Potential Salt Genes/tsnp/')) {
  snp <- read_csv(sprintf('../Potential Salt Genes/tsnp/%s', file), col_types = cols(.default = col_character()))
    if (! grepl('.', gsub(", ","",toString(colnames(snp))), fixed = TRUE)) {
      print(file)
    }
} #Only one gene reported to have no indels CHR3-11774086-11778008
#and CHR5-7474468-7477040 has no insertions



#rename
gtbl <- read_csv('../Potential Salt Genes/salty genes tr.csv', col_types = cols(.default = col_character()))
msut <-  read_csv('../Potential Salt Genes/salty msu.csv', col_types = cols(.default = col_character()))
gtbl <- merge(gtbl, msut %>% select(tr1 = msu, chr = chr, start = start, end = end, withpro = withpro), by = 'tr1')
gtbl <- gtbl[!duplicated(gtbl$tr1),]
gtbl$wGName <- unlist(lapply(gtbl$`Gene symbol synonym(s)`, str_extract, pattern = '([A-Za-z0-9])+'))
gtbl$wGName[which(is.na(gtbl$wGName))] <- paste0('NoName', c(1:length(which(is.na(gtbl$wGName)))))


new <- lapply(gtbl$tr1, function(x) {
  gname <- gtbl$wGName[gtbl$tr1 == x]
  chr <- gtbl$chr[gtbl$tr1 == x]
  start <- gtbl$withpro[gtbl$tr1 == x]
  end <- gtbl$end[gtbl$tr1 == x]
  paste(paste(chr, start, end, sep = "-"), gname, sep = ",")}
  )
gtbl$filename <- unlist(new)

old <- list.files('../Potential Salt Genes/tsnp/')
old <- lapply(list.files('../Potential Salt Genes/tsnp/'), function(x) {
  old[old == x] <- paste(strsplit(x, '-')[[1]][2], strsplit(x, '-')[[1]][3], strsplit(x, '-')[[1]][4], sep = '-')
} )

#check duplicates
files <- list.files('../Potential Salt Genes/tsnp/')
dups <- old[duplicated(old)]
old <- old[!duplicated(old)]
for (dup in dups){
  kill <- which(grepl(dup, files))[2:length(which(old == dup))]
  for (k in kill)
  {file.remove(sprintf('../Potential Salt Genes/tsnp/%s', files[k]))}
  
}


#rename 1
for (file in list.files('../Potential Salt Genes/tsnp/')) {
for (identifier in old) {
  if (grepl(identifier, file, fixed = TRUE)) {
    file.rename(from = sprintf('../Potential Salt Genes/tsnp/%s', file), to = sprintf('../Potential Salt Genes/tsnp/%s', identifier))
  }
}}

  # new <- lapply(new, function(name) { 
  # name <- gsub('\\|', '-', name)
  # name <- gsub('\\*', '(star)', name)
  # name <- gsub('\\/', '(slash)', name)
  # name <- gsub('\\:', '(colon)', name)
  # #name <- gsub('\\', '(slash2)', name)
  # } )

#rename 2
for (file in list.files('../Potential Salt Genes/tsnp/')) {
  for (name in new) {
    if (grepl(file, name[1], fixed = TRUE)) {
      file.rename(sprintf('../Potential Salt Genes/tsnp/%s', file), sprintf('../Potential Salt Genes/tsnp/%s.csv', name[1]))
    }
  }}

write_csv(gtbl, '../Potential Salt Genes/gtbl 19cols.csv')
