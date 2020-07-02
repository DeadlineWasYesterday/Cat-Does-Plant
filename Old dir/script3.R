# load the DECIPHER library in R
library(DECIPHER)
library(tidyverse)
library(seqinr)

idx_tbl <- read_csv('../Data/mf2.csv', col_types = cols(.default = col_character()))
srt_tbl <- read_csv('../Data/mf1.csv', col_types = cols(.default = col_character()))

idx_tbl2 <- read_csv('../Data/129t.csv', col_types = cols(.default = col_character()))
srt_tbl2 <- read_csv('../Data/129mf.csv', col_types = cols(.default = col_character()))

gene <- read.fasta("../Data/138genes.fasta")

sgene <- c(gene[1], gene[order(srt_tbl$SES) + 1])

sgene <- lapply(sgene, toupper)

protein <- read.fasta("../Data/129proteins.fasta")

sprotein <- c(protein[1], protein[order(srt_tbl2$SES) + 1])

sprotein <- lapply(sprotein, toupper)

gheaders <- c('original', idx_tbl$`JAPONICA NIPPONBARE POSITIONS`)
pheaders <- c('original', idx_tbl2$`JAPONICA NIPPONBARE POSITIONS`)
  
write.fasta(sgene, gheaders, ("../Data/138sgenes.fasta"), 
            open = "w", nbchar = 60000, as.string = FALSE)
write.fasta(sprotein, pheaders, ("../Data/129sproteins.fasta"), 
            open = "w", nbchar = 60000, as.string = FALSE)


# specify the path to the FASTA file (in quotes)
fas <- "../Data/129sproteins.fasta"

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readAAStringSet(fas)

# look at some of the sequences (optional)
#seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
#seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned,
                file="129saligned.fasta")



# specify the path to the FASTA file (in quotes)
fas <- "../Data/138sgenes.fasta"

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet(fas)

# look at some of the sequences (optional)
#seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned,
                file="138g_saligned.fasta")