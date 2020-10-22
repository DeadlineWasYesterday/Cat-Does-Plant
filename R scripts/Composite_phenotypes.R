library(tidyverse)

df <- read_csv('../Pie/39traits.csv')


compose <- function(traits) {
wdf <- df[c('accession_name', traits)]
tdf <- na.omit(wdf[c(traits)])
mat <- data.matrix(tdf)
pca <- prcomp(mat, scale. = TRUE)

#plot(pca$x[,1], pca$x[,2])

pca.var <- pca$sdev^2
pca.var <- round(pca.var/sum(pca.var)*100, 1)
#barplot(pca.var[1:5], main = 'Scree Plot', xlab = 'Principal Component',
#        ylab = 'Percent Variation')

pca.var[1]

}


#Seed C1
compose(c('Width_WO_husk', 'Height_WO_husk', 'Seed_Weight'))

#Seed C2
compose(c('Width_WO_husk', 'Height_WO_husk', 'Seed_Weight', 
          'Length_WO_husk'))

#Salinity C3
compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Chlorophyll_A', 
          'SES'))

#Salinity C4
compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length', 
          'Stress_Root_length',  'Stress_Chlorophyll_A', 'Stress_Chlorophyll_B',
          'SES'))

#Salinity C5
compose(c('Stress_Root_length', 'SES', 'Stress_Chlorophyll_A'))

#Salinity C6
compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Chlorophyll_A',
          'Stress_Root_length', 'SES', 'Stress_Root_Na+'))

#Salinity C7
compose(c('Stress_Chlorophyll_A', 'Stress_Root_length', 
          'SES', 'Stress_Root_Na+'))

#Biomass C8
compose(c('Shoot_weight', 'Root_weight', 'Shoot_length', 'Root_length'))

#Stress biomass C9
compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length',
          'Stress_Root_length'))

#All biomass C10
compose(c('Shoot_weight', 'Root_weight', 'Shoot_length', 'Root_length',
          'Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length',
          'Stress_Root_length'))

#Stress ions C11

compose(c('Stress_Shoot_Na+', 'Stress_Shoot_K+',
          'Stress_Root_Na+', 'Stress_Root_K+'))

#Stress Na ions C12
compose(c('Stress_Shoot_Na+', 'Stress_Root_Na+'))

#Stress Root ions C13
compose(c('Stress_Root_Na+', 'Stress_Root_K+'))

#Lost biomass C14
compose(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Shoot_length',
          'Lost_Root_length'))

#All losses C15
compose(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Shoot_length',
          'Lost_Root_length', 'Lost_Chlorophyll_A', 'Lost_Chlorophyll_B'))

#Stress induced losses C16
compose(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Chlorophyll_A',
          'SES'))

#Stress induced losses c18
compose(c('Lost_Root_length', 'Lost_Chlorophyll_A',
          'SES'))

#Biomass and chlorophyll C19
compose(c('Lost_Root_weight', 'Lost_Shoot_weight', 
          'Lost_Chlorophyll_A', 'Lost_Chlorophyll_B'))

ol <- list() 
for (i3 in 2:5) {

combs <- combn(colnames(df[2:length(df)]), i3)

odf <- data.frame()
t1 <- c()
t2 <- c()

for (i in 1:dim(combs)[2]){
  temptraits <- c()
  for (i2 in 1:dim(combs)[1]){
    temptraits <- c(temptraits, combs[i2,i])
  }
  t1 <- c(t1, paste(temptraits, collapse = ', '))
  t2 <- c(t2, compose(temptraits))
}
names(t1) <- 'traits'
names(t2) <- 'PC1 explanation'
odf <- data.frame(t1, t2)
ol[[(i3 - 1)]] <- odf

}


