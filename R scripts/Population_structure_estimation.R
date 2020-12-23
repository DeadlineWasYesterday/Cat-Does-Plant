#Population structure estimation using scatter plots, PCA and K-means clustering

list.of.packages <- c("rgl","ggplot2","knitr","rglwidget")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages,function(x){library(x,character.only=TRUE)})


library(plot3D)
library(ggpubr)
library(factoextra)
library(ggplot2)
library(extrafont)
extrafont::font_import()



pca <- read.csv('GAPIT.PCA.csv')
names(pca)[1] <- 'IRIS.code'
tab <- read.table(file = "clipboard", 
                  sep = "\t", header=TRUE)
tab['IRIS.code'] <- lapply(tab[1], str_replace, pattern = ' ', replacement = '_')


x <- merge(pca, tab, by = 'IRIS.code')

plot3d(pca$PC1, pca$PC2, pca$PC3)

plot(pca$PC1, pca$PC2, col = 'blue', xlab = 'Principal component 1', 
     ylab = 'Principal component 2', 
     main = 'Population Structure Estimation with Principal Components')

#abline(v=-50, col="purple")
#abline(h=1300, col="purple")

#abline(v=-400, col="green")
#abline(h=500, col="green")

#abline(v=-400, col="green")
#abline(h=500, col="green")

segments(x0 = -50, y0 = 1300, x1 = 300, y1 = 1300, col = 'purple', lwd = 2.5)
segments(x0 = -50, y0 = 1300, x1 = -50, y1 = 2000, col = 'purple', lwd = 2.5)
segments(x0 = 300, y0 = 1300, x1 = 300, y1 = 2000, col = 'purple', lwd = 2.5)


segments(x0 = -450, y0 = -500, x1 = -450, y1 = 500, col = 'green', lwd = 2.5)
segments(x0 = -450, y0 = 500, x1 = -1200, y1 = 500, col = 'green', lwd = 2.5)

segments(x0 = 300, y0 = -400, x1 = 300, y1 = 400, col = 'orchid1', lwd = 2.5)
segments(x0 = 300, y0 = 400, x1 = 1200, y1 = 400, col = 'orchid1', lwd = 2.5)




theme_set(theme_bw(base_size = 13, base_family = "Roboto Condensed"))

ggplot(pca, aes(x = PC1, y = PC2, label = taxa)) 
  + geom_text()
  #geom_point(color = "firebrick") +
  labs(x = "Principal component 1",
       y = "Principal component 2") + 
  coord_cartesian(xlim=c(-1100 , 1100), ylim = c(-400, 1800)) +
  geom_segment(aes(x = -50, y = 1300, xend = 300, yend = 1300), col = 'purple')  
 # geom_segment(aes(x = -50, y = 1300, xend = -50, yend = 2000), col = 'purple') + 
 # geom_segment(aes(x = 300, y = 1300, xend = 300, yend = 2000), col = 'purple') 
  



km <- kmeans( pca[2:3], centers = 3 )

fviz_cluster(km, data = pca[2:3],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_gray(), 
             main = 'Population Structure Estimation with K-Means Clustering- 2 Components',
             xlab = 'Principal component 1',
             ylab = 'Principal component 2',
             pointsize = 2.0)

fviz_cluster(km, data = pca[2:4],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_gray(), 
             main = 'Population Structure Estimation with K-Means Clustering- 3 Components',
             pointsize = 2.0)

  

##Sorting
plot(x$PC1, x$PC2, col = 'blue', xlab = 'Principal component 1', 
     ylab = 'Principal component 2', 
     main = 'Population Structure Estimation with Principal Components', )
text(x$PC1, x$PC2, labels = x$X3KRG.subpopulation)
text(x$PC1, x$PC2, labels = x$Group)
text(x$PC1, x$PC2, labels = x$IRIS.code)

x['Inferred Population'] <- ''

x[which(x$PC1 < 0 & x$PC2 < 500), 'Inferred Population'] <- 'P1'
x[which(x$PC1 > 0 & x$PC2 < 500), 'Inferred Population'] <- 'P2'
x[which(x$PC1 > -500 & x$PC2 > 1000), 'Inferred Population'] <- 'P3'

x[which(x$IRIS.code == 'IRIS_313-11223'), 'Inferred Population'] <- 'None'
x[which(x$IRIS.code == 'IRIS_313-11486'), 'Inferred Population'] <- 'None'
x[which(x$IRIS.code == 'IRIS_313-10980'), 'Inferred Population'] <- 'None'


km2 <- kmeans( x[c('PC1', 'PC2')], centers = 3 )
km3 <- kmeans( x[c('PC1', 'PC2', 'PC3')], centers = 3 )

x['Inference by K-Means with 2 components'] <- km2$cluster
x['Inference by K-Means with 3 components'] <- km3$cluster


#Collate and align (hardcoded)
fs <- read.table('WF2.3.meanQ') #fastStructure results
fs <- round(fs)

x['fastStructure Inference'] <- ''
x$`fastStructure Inference`[fs[1] == 1] <- 'P2'
x$`fastStructure Inference`[fs[2] == 1] <- 'P1'
x$`fastStructure Inference`[fs[3] == 1] <- 'P3'


x$`Inference by K-Means with 2 components`[x$`Inference by K-Means with 2 components` == 2] <- 'P1'
x$`Inference by K-Means with 2 components`[x$`Inference by K-Means with 2 components` == 1] <- 'P3'
x$`Inference by K-Means with 2 components`[x$`Inference by K-Means with 2 components` == 3] <- 'P2'

x$`Inference by K-Means with 3 components`[x$`Inference by K-Means with 3 components` == 2] <- 'P1'
x$`Inference by K-Means with 3 components`[x$`Inference by K-Means with 3 components` == 1] <- 'P3'
x$`Inference by K-Means with 3 components`[x$`Inference by K-Means with 3 components` == 3] <- 'P2'

out <- select(x, `IRIS code` = IRIS.code, `IRGC code` = IRGC.code, 
              `3KRG entry name` = X3KRG.entry.name, `3KRG subpopulation` = X3KRG.subpopulation,
              `Local grouping` = Group, `Inference by K-Means with 2 components`,
              `Inference by K-Means with 3 components`, `fastStructure Inference`,
              `Inferred population` = `Inferred Population`)
out$`Inferred population`[out$`Inferred population` == 'None'] <- 'Excluded'

write_csv(out, 'Figures and tables/S1 Table Vegetation and Population Structure.csv')
