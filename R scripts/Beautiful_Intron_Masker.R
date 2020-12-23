#given a pair of fasta files having genes and their cdses respectively,
#masks the introns from all the cdses i.e. replaces then with gaps.


library(seqinr)
library(DECIPHER)


gens <- read.fasta('Sequence and SNP data/340 candidate genes.fasta')
cds <- read.fasta('Sequence and SNP data/427 candidate cds.fasta')

noinlist <- DNAStringSet()
for (codseq in names(cds)) {
  gn <- str_extract(codseq, '[A-z0-9_]*')
  pair <- DNAStringSet(c('cds' = unlist(getSequence(cds[codseq], as.string = TRUE)),
                         'gen' = unlist(getSequence(gens[gn], as.string = TRUE))),
                       use.names = TRUE)
  pair <- AlignSeqs(pair)
  
  noint <- unname(as.character(pair['cds']))
  names(noint) <- codseq
  
  noinlist <- DNAStringSet(c(noinlist, noint))
  
}

writeXStringSet(noinlist, 'Sequence and SNP data/427extracted_genes_noint.fasta')





#recycled

#original pairer
# pair <- DNAStringSet(c('cds' = unlist(getSequence(cds[codseq], as.string = TRUE)),
#                        'gen' = unlist(getSequence(gens[str_extract(codseq, '([A-Za-z0-9_])+')], as.string = TRUE))),
#                      use.names = TRUE)