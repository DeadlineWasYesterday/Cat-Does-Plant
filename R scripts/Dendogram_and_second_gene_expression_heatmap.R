#Dendogram for clustering by gene expression and heatmap of tolerance qtls

library(ggplot2)
library(tidyverse)

et <- read_csv('Figures and tables/Expression Table rqtls added.csv')
et <- et[which(et$QTL != 'qCDP2'),]
et <- et[which(et$QTL != 'qCDP4'),]
et <- et[which(et$QTL != 'qCDP5'),]
et <- et[which(et$QTL != 'qCDP14'),]
et <- et[which(et$QTL != 'qCDP15'),]

unique(et$QTL)

wt <- et[, 4:19]
wt <- log2(wt + 1)
rownames(wt) <- et$Locus_id
wt2 <- wt[rowSums(wt) > 5,] #Used 5 for heatmap2
rhcr <- hclust(dist(wt2))
rownames(wt2)[rhcr$order]
heatmap(as.matrix(wt2),labRow = et$Locus_id, labCol = colnames(wt2))
plot(rhcr, cex=0.1, res = 1500)

png("Figures and tables/dend.png",width=6400,height=3200, res = 600)

par(cex=0.1)
plot(rhcr, labels = et$Locus_id[rhcr$order])

dev.off()

png("Figures and tables/exp3.png",width=9600,height=9600, res = 1300)

heatmap(as.matrix(wt2),labRow = rownames(wt2), labCol = colnames(wt2), cexRow = 0.05)

dev.off()


# tried out k-means
km <- kmeans(wt2, 21)

fviz_cluster(km, data = wt2,
             #palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_gray(), 
             main = 'test',
             xlab = 'Principal component 1',
             ylab = 'Principal component 2',
             pointsize = 2.0)


r <- km$cluster 
r <- et$Locus_id[order(r)]

o <- dplyr::bind_cols(r,km$cluster)

write_csv(o, 'test.csv') #k-means results


cutree(rhcr, 10)

o2 <- dplyr::bind_cols(et$Locus_id[order(cutree(rhcr, 10))], et$QTL[order(cutree(rhcr, 10))], cutree(rhcr, 10))
colnames(o2) <- c('Locus ID', 'QTL', 'Expression cluster')

write_csv(o2, 'Figures and tables/Expression dendogram 10 clusters.csv') # hierarchial clustering results written to file
