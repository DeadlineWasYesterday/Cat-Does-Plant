relate <- read_csv('../Potential Salt Genes/882 snps split4.csv')
relate$log <- -log10(relate$p)

relate <- relate[order(relate$log, decreasing = TRUE),]
relate <- select(relate, -p)

relate <- relate %>% rename('gene' = 'file', 'score' = 'log')

write.csv(relate, '../Data/882genes split at 4.csv')
