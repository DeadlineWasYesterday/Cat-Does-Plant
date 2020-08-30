library(tidyverse)

genes <- read_csv('../Potential Salt Genes/msu only.csv', col_types = cols(.default = col_character()))

for (tag in genes$`Trait Ontology`) {
  if (grepl('salt tolerance', tag, fixed = FALSE)) {
    genes[which(genes$`Trait Ontology` == tag), 'salty'] <- 'yes' }
}

#length(row.names(genes[which(genes$salty == 'yes'),]))

filtered <- genes[which(genes$salty == 'yes'),]

write_csv(filtered, '../Potential Salt Genes/salty genes.csv', na = '')
