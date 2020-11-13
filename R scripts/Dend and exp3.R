# install and load ggplot2
library(ggplot2)
library(tidyverse)

# import the data
ld <- read.table("../plink.ld", sep="\t",header=T)

# plot the average correlation for each snp distance
ggplot(ld) +
  geom_line(aes(x=BP_B - BP_A, y = R2))


x <- read_csv('Notebooks/Exp1.csv')
x2 <- x[221:477, 1:17]
rhcr <- hclust(dist(x2[2:17]))
x2$Locus_id[rhcr$order]
heatmap(as.matrix(x2[2:17]),labRow = x2$Locus_id, labCol = colnames(x2))

x <- read_csv('Notebooks/Exp2.csv')
x2 <- x[221:477, 1:33]
rhcr <- hclust(dist(x2[2:33]))
x2$Locus_id[rhcr$order]
heatmap(as.matrix(x2[2:33]),labRow = x2$Locus_id, labCol = colnames(x2[2:33]), mar=c(14,14))


x3 <- x2[rhcr$order,][115:257, 1:33]
rhcr <- hclust(dist(x3[2:33]))
x3$Locus_id[rhcr$order]
heatmap(as.matrix(x3[2:33]),labRow = x3$Locus_id, labCol = colnames(x3[2:33]), mar=c(14,14))


et <- read_csv('Figures and tables/Expression Table rqtls added.csv')
# et <- et[which(et$QTL != 'qCDP2'),]
# et <- et[which(et$QTL != 'qCDP4'),]
# et <- et[which(et$QTL != 'qCDP5'),]
# et <- et[which(et$QTL != 'qCDP14'),]
# et <- et[which(et$QTL != 'qCDP15'),]

unique(et$QTL)

wt <- et[, 4:19]
wt <- log2(wt + 1)
rownames(wt) <- et$Locus_id
wt2 <- wt[rowSums(wt) > 1,] #Used 5 for heatmap
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



km <- kmeans(wt2, 21)

fviz_cluster(km, data = wt2,
             #palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_gray(), 
             main = 'Population Structure Estimation with K-Means Clustering- 2 Components',
             xlab = 'Principal component 1',
             ylab = 'Principal component 2',
             pointsize = 2.0)


r <- km$cluster 
r <- et$Locus_id[order(r)]

o <- dplyr::bind_cols(r,km$cluster)

write_csv(o, 'test.csv')


cutree(rhcr, 10)

o2 <- dplyr::bind_cols(et$Locus_id[order(cutree(rhcr, 10))], et$QTL[order(cutree(rhcr, 10))], cutree(rhcr, 10))
colnames(o2) <- c('Locus ID', 'QTL', 'Expression cluster')

write_csv(o2, 'Figures and tables/Expression dendogram 10 clusters.csv')
