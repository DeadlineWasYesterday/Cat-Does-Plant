library(tidyverse)

#Make sure all tables have snps
for (file in list.files('../Potential Salt Genes/snpia/')) {
  snp <- read_csv(sprintf('../Potential Salt Genes/snpia/%s', file), col_types = cols(.default = col_character()))
    if (! grepl('.', gsub(", ","",toString(colnames(snp))), fixed = TRUE)) {
      print(file)
    }
} #Only one gene reported to have no indels CHR3-11774086-11778008
#and CHR5-7474468-7477040 has no insertions

#rename
gtbl <- read_csv('../Potential Salt Genes/salty genes tr.csv', col_types = cols(.default = col_character()))

msut <-  read_csv('../Potential Salt Genes/salty msu.csv', col_types = cols(.default = col_character()))

gtbl <- merge(gtbl, msut %>% rename(tr1 = msu), by = 'tr1')

new <- lapply(gtbl$tr1, function(x) {
  gname <- gtbl$wGName[gtbl$tr1 == x]
  chr <- gtbl$chr[gtbl$tr1 == x]
  start <- gtbl$withpro[gtbl$tr1 == x]
  end <- gtbl$end[gtbl$tr1 == x]
  paste(chr, start, end, gname, sep = "-")}
  )
gtbl$filename <- new

old <- list.files('../Potential Salt Genes/snpia/')
old <- lapply(list.files('../Potential Salt Genes/snpia/'), function(x) {
  old[old == x] <- paste(strsplit(x, '-')[[1]][2], strsplit(x, '-')[[1]][3], strsplit(x, '-')[[1]][4], sep = '-')
} )

#check duplicates
old[duplicated(old)]

#rename 1
for (file in list.files('../Potential Salt Genes/snpia/')) {
for (identifier in old) {
  if (grepl(identifier, file, fixed = TRUE)) {
    file.rename(from = sprintf('../Potential Salt Genes/snpia/%s', file), to = sprintf('../Potential Salt Genes/snpia/%s', identifier))
  }
}}

new <- lapply(new, function(name) { 
  name <- gsub('\\|', '-', name)
  name <- gsub('\\*', '(star)', name)
  name <- gsub('\\/', '(slash)', name)
  name <- gsub('\\:', '(colon)', name)
  #name <- gsub('\\', '(slash2)', name)
  } )

#rename 2
for (file in list.files('../Potential Salt Genes/snpia/')) {
  for (name in new) {
    if (grepl(file, name[1], fixed = TRUE)) {
      file.rename(sprintf('../Potential Salt Genes/snpia/%s', file), sprintf('../Potential Salt Genes/snpia/%s.csv', name[1]))
    }
  }}

write_csv(gtbl[1:19], '../Potential Salt Genes/gtbl 19cols.csv')
