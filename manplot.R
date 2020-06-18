# Load the library
library(qqman)

results <- read_csv('../Potential Salt Genes/882 snps split8.csv')
relate <- read_csv('../Potential Salt Genes/882 snps split8.csv')
relate$log <- -log10(relate$p)

results$file <- as.integer(unlist(lapply(results$file, str_extract, pattern = "[[0-9]]+")))
results <- results %>% rename('chr' = 'file')

loc = c()
for (cn in unique(results$chr)) {
  loc = c(loc, 1:sum(results$chr == cn))
}

results$loc <- loc
results <- results %>% rename('snp' = 'location')

# Make the Manhattan plot on the gwasResults dataset
manhattan(results, chr="chr", bp="loc", snp="snp", p="p",
          col = c('#A93226', '#1F618D', '#1E8449', '#283747'),
          #col = c('#9B59B6', '#2E86C1', '#17A589', '#28B463'),
          genomewideline = 4,
          annotateTop = True)

