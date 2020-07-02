library(tidyverse)

df <- read_csv('../Potential Salt Genes/bfch2.csv', col_types = cols(.default = col_character()))

df <- df[c(1,5,31:length(df))]

for (col in 3:length(df)) {
  #fill NA with Ref
  df[col][is.na(df[col])] <- unlist(df[col])[1]
  #change anything that isn't ref into 1
  df[col][df[col] != unlist(df[col])[1]] <- '1'
  #change ref to 0
  df[col][df[col] == unlist(df[col])[1]] <- '0'
  #change to numneric
  df[col] <- as.numeric(unlist(df[col]))
  
}


mat <- data.matrix(df[3:length(df)])
rownames(mat) <- c('original', unname(unlist(na.omit(df[2]))))
plants <- unname(unlist(na.omit(df[1])))
subpops <- c('original', unname(unlist(na.omit(df[2]))))

pca <- prcomp(mat, scale. = TRUE)

plot(pca$x[,1], pca$x[,2])

pca.var <- pca$sdev^2
pca.var <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var[1:10], main = 'Scree Plot', xlab = 'Principal Component',
        ylab = 'Percent Variation')


library(ggplot2)
pca.data <- data.frame(Sample=rownames(mat),
                       X=pca$x[,1],
                       Y=pca$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=subpops)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var[2], "%", sep="")) +
  theme_bw() +
  ggtitle("2183Genes CHR2 PC1 vs PC2")

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)