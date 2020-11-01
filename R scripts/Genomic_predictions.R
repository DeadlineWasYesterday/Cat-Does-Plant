library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl<-makeCluster(30, type = "SOCK")
registerDoSNOW(cl)
library(rrBLUP)
library(dplyr)

myY <- read.table('../../Data/39transformed_WF2.txt', header = TRUE)
myG <- read.table('../../GWAS_factory/GnQWF2/GAPIT.Genotype.Numerical.txt', header = TRUE)

M <- as.matrix(myG[,2:length(myG)])

myO <- myY[,1]

for (i in 2:26) {
  blist <- c()
  for (i2 in 1:1) {
  pred <- mixed.solve(y = myY[,i], M)
  
  blist <- c(blist, pred$beta)
  }
  BLUE <- mean(blist)
  
  r <- M %*% pred$u
  r <- r[,1] + BLUE
  names(r) <- names(myY[,i])
  
  myO <- dplyr::bind_cols(myO, r)
  
}


# #random population of 200 lines with 1000 markers
# M <- matrix(rep(0,200*1000),200,1000)
# for (i in 1:200) {
#   M[i,] <- ifelse(runif(1000)<0.5,-1,1)
# }
# #random phenotypes
# u <- rnorm(1000)
# g <- as.vector(crossprod(t(M),u))
# h2 <- 0.5 #heritability
# y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
# #predict marker effects
# ans <- mixed.solve(y,Z=M) #By default K = I
# accuracy <- cor(u,ans$u)
# #predict breeding values
# ans <- mixed.solve(y,K=A.mat(M))
# accuracy <- cor(g,ans$u)
# cor(y, ans$u)
# z