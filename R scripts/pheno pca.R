library(tidyverse)

df <- read_csv('../Data/GWAS data/163x23 both.csv')

df <- df[c(6,8,10,12,20,22,23)]

df <- scale(df)
mat <- data.matrix(df)

pca <- prcomp(mat, scale. = TRUE)

plot(pca$x[,1], pca$x[,2])


pca.var <- pca$sdev^2
pca.var <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var[1:10], main = 'Scree Plot', xlab = 'Principal Component',
        ylab = 'Percent Variation')

write_csv(df, '../Data/GWAS data/163x25 wpca.csv', na = '')


library(tidyverse)

df <- read_csv('../Data/GWAS data/183x32 both.csv')

df <- df[c(7,9,11,13,25,27,32)]

df <- unname(df)

idf <- mice(df, m=5, maxit = 50, method = 'pmm')
cdf <- complete(idf, 2)

#df <- scale(df)
mat <- data.matrix(cdf)

pca <- prcomp(mat, scale. = TRUE)

plot(pca$x[,1], pca$x[,2])


pca.var <- pca$sdev^2
pca.var <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var[1:10], main = 'Scree Plot', xlab = 'Principal Component',
        ylab = 'Percent Variation')

odf <- read_csv('../Data/GWAS data/183x32 both.csv')
odf <- dplyr::bind_cols(odf,pca$x[,1],pca$x[,2])

write_csv(odf, '../Data/GWAS data/183x34 ipc.csv', na = '')


#No imputation

library(tidyverse)

df <- read_csv('../Data/GWAS data/183x32 both.csv')
df <- df[c(1,7,9,11,13,25,27,32)]
df <- na.omit(df)

#df <- scale(df)
mat <- data.matrix(df)

pca <- prcomp(mat, scale. = TRUE)

plot(pca$x[,1], pca$x[,2])


pca.var <- pca$sdev^2
pca.var <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var[1:10], main = 'Scree Plot', xlab = 'Principal Component',
        ylab = 'Percent Variation')


df <- dplyr::bind_cols(df,pca$x[,1],pca$x[,2])

odf <- read_csv('../Data/GWAS data/183x32 both.csv')
odf <- merge(odf,df[c(1,9,10)])

write_csv(odf, '../Data/GWAS data/173x34 pc.csv', na = '')
