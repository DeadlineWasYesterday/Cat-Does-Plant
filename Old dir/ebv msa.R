library(DECIPHER)

# specify the path to the FASTA file (in quotes)
fas <- "../Data/ebv for msa"

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet(fas)

# look at some of the sequences (optional)
#seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
#seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=0)
