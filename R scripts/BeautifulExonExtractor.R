library(seqinr)
library(DECIPHER)


gens <- read.fasta('../Data/RAP-DB/Putative RAPDB genes.fasta')
cds <- read.fasta('../Data/RAP-DB/Putative RAPDB cds.fasta')

exlist <- DNAStringSet()
for (codseq in names(cds)) {
  pair <- DNAStringSet(c('cds' = unlist(getSequence(cds[codseq], as.string = TRUE)),
                'gen' = unlist(getSequence(gens[codseq], as.string = TRUE))),
               use.names = TRUE)
  pair <- AlignSeqs(pair)
  
  exs <- unlist(str_extract_all(pair$cds, '([A-Z])+'))
  names(exs) <- paste0(paste0(codseq, '_exon_'), (1:length(exs)))
  
  exlist <- DNAStringSet(c(exlist, exs))
  
}

writeXStringSet(exlist, '../Data/RAP-DB/RAPDB Putative exons.fasta')





#recycled

#original pairer
# pair <- DNAStringSet(c('cds' = unlist(getSequence(cds[codseq], as.string = TRUE)),
#                        'gen' = unlist(getSequence(gens[str_extract(codseq, '([A-Za-z0-9_])+')], as.string = TRUE))),
#                      use.names = TRUE)