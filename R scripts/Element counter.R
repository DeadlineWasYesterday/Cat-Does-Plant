#Count elements in upstream gene regions of Oryza

library(stringr)
library(seqinr)

#Load 1k, 2k and 3k upstream gene sequences
dna1 <- read.fasta('../Data/RAP-DB/IRGSP-1.0_1kb-upstream_2020-06-03.fasta')

dna2 <- read.fasta('../Data/RAP-DB/IRGSP-1.0_2kb-upstream_2020-06-03.fasta')

dna3 <- read.fasta('../Data/RAP-DB/IRGSP-1.0_3kb-upstream_2020-06-03.fasta')


#Count GATACA in 1k upstream promoters
gataca.in1kup <- str_count(sapply(dna1, c2s), "gataca")
names(gataca.in1kup) <- names(dna1)

#Count ttgaca in 2k upstream promoters
ttgaca.in2kup <- str_count(sapply(dna2, c2s), "ttgaca")
names(ttgaca.in2kup) <- names(dna2)

#Count ggccaatct in 3k upstream promoters
ggccaatct.in3kup <- str_count(sapply(dna3, c2s), "ggccaatct")
names(ggccaatct.in3kup) <- names(dna3)


#save result

table <- cbind(data.frame(gataca.in1kup), data.frame(ttgaca.in2kup), 
               data.frame(ggccaatct.in3kup))

write.csv(table,'../Data/Element Counts.csv', na = '')

