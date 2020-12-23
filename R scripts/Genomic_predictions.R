#genomic predictions using mixed.solve

library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl<-makeCluster(30, type = "SOCK")
registerDoSNOW(cl)
library(rrBLUP)
library(dplyr)
library(tidyverse)

#predictions of observed phenotypes
myY <- read_csv('39traitsord.csv')
myG <- read_csv('../WF2m05.mat.csv', col_names = FALSE)
M <- t(myG)

myO <- as.data.frame(myY[,1]) 
names(myO) <- 'accession_name'

for (i in 2:40) {
  M <- t(myG)[!is.na(myY[,i]),]
  
  pred <- mixed.solve(y = as.matrix(na.omit(myY[,i])), as.matrix(M))

  BLUE <- pred$beta
  
  r <- M %*% pred$u
  r <- as.vector(r) + as.numeric(BLUE)
  
  myA <- dplyr::bind_cols(myY[,1][!is.na(myY[,i])], r)
  names(myA) <- names(myY)[c(1,i)]
  myO <- merge(myO, myA, by = 'accession_name', all = TRUE)
  
  print(head(myO))
  
}

write_csv(myO, 'PredSNPeffects raw.csv')


#predictions after transformation by warpedlmm
myY <- read_csv('39transbyWF2ord.txt')
myG <- read_csv('../WF2m05.mat.csv', col_names = FALSE)
M <- t(myG)

myO <- as.data.frame(myY[,1]) 
names(myO) <- 'accession_name'

for (i in 2:40) {
  M <- t(myG)[!is.na(myY[,i]),]
  
  pred <- mixed.solve(y = as.matrix(na.omit(myY[,i])), as.matrix(M))

  BLUE <- pred$beta
  
  r <- M %*% pred$u
  r <- as.vector(r) + as.numeric(BLUE)
  
  myA <- dplyr::bind_cols(myY[,1][!is.na(myY[,i])], r)
  names(myA) <- names(myY)[c(1,i)]
  myO <- merge(myO, myA, by = 'accession_name', all = TRUE)
  
  print(head(myO))
  
}

write_csv(myO, 'PredSNPeffects transformed.csv')