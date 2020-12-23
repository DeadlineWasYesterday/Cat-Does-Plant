#a deprecated file for making composite phenotypes by 
#multidimensional scaling. 

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



a <- c()
b <- c()

#Seed C1
a <- append(a, paste(c('Width_WO_husk', 'Height_WO_husk', 'Seed_Weight'), collapse = ', '))
b <- append(b, compose(c('Width_WO_husk', 'Height_WO_husk', 'Seed_Weight')))

#Seed C2
compose(c('Width_WO_husk', 'Height_WO_husk', 'Seed_Weight', 
          'Length_WO_husk'))

a <- append(a, paste(c('Width_WO_husk', 'Height_WO_husk', 'Seed_Weight', 
                  'Length_WO_husk'), collapse = ', '))
b <- append(b, compose(c('Width_WO_husk', 'Height_WO_husk', 'Seed_Weight', 
                    'Length_WO_husk')))

#Salinity C3
compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Chlorophyll_A', 
          'SES'))

a <- append(a, paste(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Chlorophyll_A', 
                  'SES'), collapse = ', '))
b <- append(b, compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Chlorophyll_A', 
                    'SES')))


#Salinity C4
compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length', 
          'Stress_Root_length',  'Stress_Chlorophyll_A', 'Stress_Chlorophyll_B',
          'SES'))

a <- append(a, paste(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length', 
                  'Stress_Root_length',  'Stress_Chlorophyll_A', 'Stress_Chlorophyll_B',
                  'SES'), collapse = ', '))
b <- append(b, compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length', 
                    'Stress_Root_length',  'Stress_Chlorophyll_A', 'Stress_Chlorophyll_B',
                    'SES')))

#Salinity C5
compose(c('Stress_Root_length', 'SES', 'Stress_Chlorophyll_A'))

a <- append(a, paste(c('Stress_Root_length', 'SES', 'Stress_Chlorophyll_A'), collapse = ', '))
b <- append(b, compose(c('Stress_Root_length', 'SES', 'Stress_Chlorophyll_A')))

#Salinity C6
compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Chlorophyll_A',
          'Stress_Root_length', 'SES', 'Stress_Root_Na+'))

a <- append(a, paste(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Chlorophyll_A',
                  'Stress_Root_length', 'SES', 'Stress_Root_Na+'), collapse = ', '))
b <- append(b, compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Chlorophyll_A',
                    'Stress_Root_length', 'SES', 'Stress_Root_Na+')))

#Salinity C7
compose(c('Stress_Chlorophyll_A', 'Stress_Root_length', 
          'SES', 'Stress_Root_Na+'))

a <- append(a, paste(c('Stress_Chlorophyll_A', 'Stress_Root_length', 
                  'SES', 'Stress_Root_Na+'), collapse = ', '))
b <- append(b, compose(c('Stress_Chlorophyll_A', 'Stress_Root_length', 
                    'SES', 'Stress_Root_Na+')))

#Biomass C8
compose(c('Shoot_weight', 'Root_weight', 'Shoot_length', 'Root_length'))

a <- append(a, paste(c('Shoot_weight', 'Root_weight', 'Shoot_length', 'Root_length'), collapse = ', '))
b <- append(b, compose(c('Shoot_weight', 'Root_weight', 'Shoot_length', 'Root_length')))

#Stress biomass C9
compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length',
          'Stress_Root_length'))

a <- append(a, paste(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length',
                  'Stress_Root_length'), collapse = ', '))
b <- append(b, compose(c('Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length',
                    'Stress_Root_length')))

#All biomass C10
compose(c('Shoot_weight', 'Root_weight', 'Shoot_length', 'Root_length',
          'Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length',
          'Stress_Root_length'))

a <- append(a, paste(c('Shoot_weight', 'Root_weight', 'Shoot_length', 'Root_length',
                  'Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length',
                  'Stress_Root_length'), collapse = ', '))
b <- append(b, compose(c('Shoot_weight', 'Root_weight', 'Shoot_length', 'Root_length',
                    'Stress_Shoot_weight', 'Stress_Root_weight', 'Stress_Shoot_length',
                    'Stress_Root_length')))

#Stress ions C11

compose(c('Stress_Shoot_Na+', 'Stress_Shoot_K+',
          'Stress_Root_Na+', 'Stress_Root_K+'))

a <- append(a, paste(c('Stress_Shoot_Na+', 'Stress_Shoot_K+',
                  'Stress_Root_Na+', 'Stress_Root_K+'), collapse = ', '))
b <- append(b, compose(c('Stress_Shoot_Na+', 'Stress_Shoot_K+',
                    'Stress_Root_Na+', 'Stress_Root_K+')))


#Stress Na ions C12
compose(c('Stress_Shoot_Na+', 'Stress_Root_Na+'))

a <- append(a, paste(c('Stress_Shoot_Na+', 'Stress_Root_Na+'), collapse = ', '))
b <- append(b, compose(c('Stress_Shoot_Na+', 'Stress_Root_Na+')))

#Stress Root ions C13
compose(c('Stress_Root_Na+', 'Stress_Root_K+'))

a <- append(a, paste(c('Stress_Root_Na+', 'Stress_Root_K+'), collapse = ', '))
b <- append(b, compose(c('Stress_Root_Na+', 'Stress_Root_K+')))

#Lost biomass C14
compose(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Shoot_length',
          'Lost_Root_length'))

a <- append(a, paste(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Shoot_length',
                  'Lost_Root_length'), collapse = ', '))
b <- append(b, compose(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Shoot_length',
                    'Lost_Root_length')))

#All losses C15
compose(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Shoot_length',
          'Lost_Root_length', 'Lost_Chlorophyll_A', 'Lost_Chlorophyll_B'))

a <- append(a, paste(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Shoot_length',
                  'Lost_Root_length', 'Lost_Chlorophyll_A', 'Lost_Chlorophyll_B'), collapse = ', '))
b <- append(b, compose(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Shoot_length',
                    'Lost_Root_length', 'Lost_Chlorophyll_A', 'Lost_Chlorophyll_B')))

#Stress induced losses C16
compose(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Chlorophyll_A',
          'SES'))

a <- append(a, paste(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Chlorophyll_A',
                  'SES'), collapse = ', '))
b <- append(b, compose(c('Lost_Shoot_weight', 'Lost_Root_weight', 'Lost_Chlorophyll_A',
                    'SES')))

#Stress induced losses c17
compose(c('Lost_Root_length', 'Lost_Chlorophyll_A',
          'SES'))

a <- append(a, paste(c('Lost_Root_length', 'Lost_Chlorophyll_A',
                  'SES'), collapse = ', '))
b <- append(b, compose(c('Lost_Root_length', 'Lost_Chlorophyll_A',
                    'SES')))

#Biomass and chlorophyll C18
compose(c('Lost_Root_weight', 'Lost_Shoot_weight', 
          'Lost_Chlorophyll_A', 'Lost_Chlorophyll_B'))

a <- append(a, paste(c('Lost_Root_weight', 'Lost_Shoot_weight', 
                  'Lost_Chlorophyll_A', 'Lost_Chlorophyll_B'), collapse = ', '))
b <- append(b, compose(c('Lost_Root_weight', 'Lost_Shoot_weight', 
                    'Lost_Chlorophyll_A', 'Lost_Chlorophyll_B')))


odt <- data.frame(a, b)
colnames(odt) <- c('t1', 't2')



ol <- list() 
for (i3 in 2:6) {

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

o2 <- read.csv('2.csv') 
o3 <- read.csv('3.csv') 
o4 <- read.csv('4.csv') 
o5 <- read.csv('5.csv') 
o6 <- read.csv('6.csv') 

o2 <- o2[order(o2$t2, decreasing = TRUE),]
o3 <- o3[order(o3$t2, decreasing = TRUE),]
o4 <- o4[order(o4$t2, decreasing = TRUE),]
o5 <- o5[order(o5$t2, decreasing = TRUE),]
o6 <- o6[order(o6$t2, decreasing = TRUE),]

odt <- dplyr::bind_rows(odt, o2[1:21,], o2[1:10,], o4[1:10,], o5[1:10,], o6[1:10,])

colnames(odt) <- c('Trait combination', 'Explained variance')

odt <- odt[!duplicated(odt$`Trait combination`),]

write_csv(odt, 'Composite_trait_combinations.csv')





forge <- function(traits) {
  wdf <- df[c('accession_name', traits)]
  adf <- na.omit(wdf)
  tdf <- na.omit(wdf[c(traits)])
  mat <- data.matrix(tdf)
  pca <- prcomp(mat, scale. = TRUE)
  
  #plot(pca$x[,1], pca$x[,2])
  
  pca.var <- pca$sdev^2
  pca.var <- round(pca.var/sum(pca.var)*100, 1)
  #barplot(pca.var[1:5], main = 'Scree Plot', xlab = 'Principal Component',
  #        ylab = 'Percent Variation')
  
  adf <- dplyr::bind_cols(adf$accession_name, pca$x[,1])
  
  adf
  
}


cdf <- forge(strsplit(odt$`Trait combination`, ', ')[[1]])
names(cdf) <- c('accession_name', sprintf('C%d', 1))

for (i in 2:68) {
  x <- forge(strsplit(odt$`Trait combination`, ', ')[[i]])
  names(x) <- c('accession_name', sprintf('C%d', i))
  
  cdf <- merge(cdf, x, by = 'accession_name', all = TRUE)
  
}


write_csv(cdf, '68composite.csv')
write.table(cdf, '68composite.txt', row.names = FALSE, sep = '\t')

