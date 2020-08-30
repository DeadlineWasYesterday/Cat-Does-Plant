library(tidyverse)
library(ggplot2)

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
  df[col] <- as.numeric(unlist(df[col])) }


mat <- data.matrix(df[3:length(df)])
rownames(mat) <- c('original', unname(unlist(na.omit(df[2]))))
plants <- unname(unlist(na.omit(df[1])))
subpops <- c('original', unname(unlist(na.omit(df[2]))))

pca <- prcomp(mat, scale. = TRUE)

pca.var <- pca$sdev^2
pca.var <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var[1:10], main = 'Scree Plot', xlab = 'Principal Component',
        ylab = 'Percent Variation')

pca.data <- data.frame(Sample=rownames(mat),
                       X=pca$x[,1],
                       Y=pca$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=subpops)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var[2], "%", sep="")) +
  theme_bw() +
  ggtitle("2183Genes CHR2 PC1 vs PC2")



gmat <- t(unname(mat))

gmat <- sweep(gmat, 1, rowMeans(gmat))

dev <- sqrt((1 + rowSums(gmat)) / (2 + 2 * ncol(gmat)) * (1 - (1 + rowSums(gmat)) / (2 + 2 * ncol(gmat))))

gmat <- sweep(gmat, 1, dev, "/")

cmat <- cov(gmat, gmat)

emat <- eigen(cmat)


ggplot(data=pca.data, aes(x=emat$vectors[,1], y=emat$vectors[,2], label=subpops)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var[2], "%", sep="")) +
  theme_bw() +
  ggtitle("2183Genes CHR2 PC1 vs PC2")




###############################
Gmatrix.me <- Gmatrix(mat, method="Yang", maf=0.05)

eig <- eigen(Gmatrix.me)

ggplot(data=pca.data, aes(x=eig$vectors[,1], y=eig$vectors[,2], label=subpops)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var[2], "%", sep="")) +
  theme_bw() +
  ggtitle("2183Genes CHR2 PC1 vs PC2")


############################################
library(AGHmatrix)
G <- Gmatrix(as.matrix(mat), 'Yang')

mat <- data.matrix(df[3:length(df)])
rownames(mat) <- c('original', unname(unlist(na.omit(df[2]))))
plants <- unname(unlist(na.omit(df[1])))
subpops <- c('original', unname(unlist(na.omit(df[2]))))

x <- matrix(ncol = len)
for (i in 1:5) {
  step1 <- test[i,] - (sum(test[i,]) / 4)
  pi <- (1+ (sum(test[i,]))) / 10
  step2 <- step1 / sqrt(pi*(1 - pi))
  x <- rbind(x, step2)
}

x <- x[-1,]  

t2 <- matrix(ncol = 4)
for (i in 1:4) {
  t1 <- c()
  for (i2 in 1:4) {
    s <- cov(x[,i], x[,i2])
    t1 <- c(t1,s)
  }
  t2 <- rbind(t2, t1)
}

t2 <- t2[-1,]
X <- eigen(t2)


test <- matrix(c(0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0), 6,4, byrow = TRUE)
test <- test[-1,]

test2 <- sweep(test, 1, rowMeans(test))
dev <- sqrt((1 + rowSums(test)) / (2 + 2 * ncol(test)) * (1 - (1 + rowSums(test)) / (2 + 2 * ncol(test))))
test2 <- sweep(test2, 1, dev, "/")


testt <- t(test)
pca <- prcomp(testt)

pca.var <- pca$sdev^2
pca.var <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var[1:10], main = 'Scree Plot', xlab = 'Principal Component',
        ylab = 'Percent Variation')


x <- matrix(ncol = 4)
for (i in 1:6) {
  step1 <- test[i,] - (sum(test[i,]) / 4)
  pi <- (1+ (sum(test[i,]))) / 10
  step2 <- step1 / sqrt(pi*(1 - pi))
  x <- rbind(x, step2)
}

x <- x[-1,]  

t2 <- matrix(ncol = 4)
for (i in 1:4) {
  t1 <- c()
  for (i2 in 1:4) {
    s <- cov(x[,i], x[,i2])
    t1 <- c(t1,s)
  }
  t2 <- rbind(t2, t1)
}

t2 <- t2[-1,]
X <- eigen(t2)


library(ggplot2)
pca.data <- data.frame(X=pca$x[,1],
                       Y=pca$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=c(1,2,3,4))) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var[2], "%", sep="")) +
  theme_bw() +
  ggtitle("2183Genes CHR2 PC1 vs PC2")

pca.data <- data.frame(X=X$vectors[,1],
                       Y=X$vectors[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=c(1,2,3,4))) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var[2], "%", sep="")) +
  theme_bw() +
  ggtitle("2183Genes CHR2 PC1 vs PC2")
