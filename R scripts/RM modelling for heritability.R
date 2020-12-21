library(lme4)
library(tidyverse)
library(jtools)
library(lmerTest)

df <- read_csv('Notebooks/12traits for rem.csv')

df['accession_name'] <- as.factor(pull(df['accession_name']))
df['Treatment'] <- as.factor(pull(df['Treatment']))

an <- pull(df['accession_name'])
tm <- pull(df['Treatment'])

hv <- c()
for (i in 3:14) {
  tv <- pull(df[,i])
  tdf <- dplyr::bind_cols(an,tm,tv)
  colnames(tdf) <- c('an', 'tm', 'tv')
 model <- lmer(tv~(1|an)+(1|tm)+(1|an:tm),tdf) 
 tr <- as.data.frame(VarCorr(model,comp='vcov'))
 print(colnames(df)[i])
 print(tr)
 h2 <- tr[2,4] / (tr[2,4] + (tr[4,4] / 3))
 #h2 <- tr[1,4] / (tr[1,4] + (tr[3,4] / 3))
 
 hv <-c(hv, h2)
 
}

ht <- dplyr::bind_cols(colnames(df)[3:14], hv)
colnames(ht) <- c('Trait', 'Heritability')
write_csv(ht, 'Figures and tables/Tables/RE heritability.csv')

coef(model)
summ(model)
ranova(model)
