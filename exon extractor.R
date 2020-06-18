library(tidyverse)
library(seqinr)
library(DECIPHER)

msut <-  read_csv('../Potential Salt Genes/salty msu.csv', col_types = cols(.default = col_character()))

genes <- read.fasta('../Pie/884 salty genes.fasta', seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE, whole.header = TRUE)
cdses <- read.fasta('../Pie/884 salty cdses.fasta', seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE, whole.header = TRUE)


skeleton <- data.frame()
for (i in 1:length(cdses)) {
  write.fasta(c(genes[i], cdses[i]), 
              c(getName(genes[i]), getName(cdses[i])), 
              nbchar = getLength(genes[i]),
              file.out = 'temp.fasta', open = 'w')
  
  aligned <- AlignSeqs(readDNAStringSet('temp.fasta'))
  writeXStringSet(aligned, filepath = 'temp.fasta',)
  
  wsq <- read.fasta('temp.fasta', seqtype = "DNA", forceDNAtolower = FALSE, whole.header = TRUE)
  wsq <- wsq[[2]]
  
  exs <- c()
  exe <- c()
  
  pin <- 0
  while (sum(which(wsq != '-')) != 0) {
  start <- which(wsq != '-')[1]
  wsq <- wsq[start:length(wsq)]
  exs <- c(exs, start + pin)
  pin <- start + pin - 1
  end <- which(wsq == '-')[1]
  if (sum(which(wsq == '-')) == 0) { end <- length(wsq) 
  exe <- c(exe, end + pin) 
  break }
  wsq <- wsq[end:length(wsq)]
  exe <- c(exe, end + pin - 1)
  pin <- end + pin - 1
  }
  
mat <- c()
columns <- c()
for (num in 1:length(exs)) {
  mat <- c(mat, exs[num], exe[num])
  columns <- c(columns, sprintf("ex%d s", num), sprintf("ex%d e", num))
}
bone <- data.frame(t(mat))
colnames(bone) <- columns

  skeleton <- dplyr::bind_rows(skeleton, bone)
}

indices <- lapply(getName(genes), function(x)
  x <- stringr::str_remove(x, " genomic sequence "))

indices[716] <- 'LOC_Os07g33954 | OsGPT2 | Glc-6-phosphate translocator2 duped'

rownames(skeleton) <- indices

write.csv(skeleton, '../Potential Salt Genes/exon list.csv')
