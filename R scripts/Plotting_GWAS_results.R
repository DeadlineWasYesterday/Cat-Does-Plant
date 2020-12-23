#Using qqman to plot GWAS results

library(qqman)
library(tidyverse)

makeplots <- function(sel) { 

for (file in sel) {
  df <- read_csv(sprintf('../../../GWAS_factory/WF2noFDR/%s', file))
  
  name <- paste0(unlist(str_split(file, "\\."))[2:3], collapse = '.')
  
  jpeg(sprintf('%s.jpg', name), width = 1366, height = 788)
  
  manhattan(df, chr="Chromosome", bp="Position", snp="SNP", p="P.value",
            col = c('#9B59B6', '#2E86C1', '#17A589', '#28B463'),
            genomewideline = 6.5,
            annotateTop = True)
  
  dev.off()
}
}

fv <- str_extract(dir('../../../GWAS_factory/WF2noFDR/'), '.*\\Blink.*GWAS.Results.csv')
fv <- fv[!is.na(fv)]

makeplots(fv)

WF2sig <- c('GAPIT.Blink.Chlorophyll_B.GWAS.Results.csv',
            'GAPIT.Blink.Height_WO_husk.GWAS.Results.csv',
            'GAPIT.Blink.Length_WO_husk.GWAS.Results.csv',
            'GAPIT.Blink.Lost_Shoot_length.GWAS.Results.csv',
            'GAPIT.Blink.Seed_volume.GWAS.Results.csv',
            'GAPIT.Blink.Shoot_K..GWAS.Results.csv',
            'GAPIT.CMLM.Chlorophyll_B.GWAS.Results.csv',
            'GAPIT.CMLM.Lost_Shoot_length.GWAS.Results.csv',
            'GAPIT.CMLM.Shoot_K..GWAS.Results.csv')

makeplots(WF2sig)

WF2sug <- c('GAPIT.Blink.Chlorophyll_A.GWAS.Results.csv', 
            'GAPIT.Blink.Chlorophyll_B.GWAS.Results.csv', 
            'GAPIT.Blink.Height_WO_husk.GWAS.Results.csv', 
            'GAPIT.Blink.Length_WO_husk.GWAS.Results.csv', 
            'GAPIT.Blink.Lost_Chlorophyll_B.GWAS.Results.csv', 
            'GAPIT.Blink.Lost_Root_length.GWAS.Results.csv', 
            'GAPIT.Blink.Lost_Root_thickness.GWAS.Results.csv', 
            'GAPIT.Blink.Lost_Root_weight.GWAS.Results.csv', 
            'GAPIT.Blink.Lost_Shoot_length.GWAS.Results.csv', 
            'GAPIT.Blink.Lost_Shoot_thickness.GWAS.Results.csv', 
            'GAPIT.Blink.Root_K..GWAS.Results.csv', 
            'GAPIT.Blink.Root_thickness.GWAS.Results.csv', 
            'GAPIT.Blink.Seed_density.GWAS.Results.csv', 
            'GAPIT.Blink.Seed_volume.GWAS.Results.csv', 
            'GAPIT.Blink.Seed_Weight.GWAS.Results.csv', 
            'GAPIT.Blink.Shoot_K..GWAS.Results.csv', 
            'GAPIT.Blink.Shoot_Na..GWAS.Results.csv', 
            'GAPIT.Blink.Shoot_thickness.GWAS.Results.csv', 
            'GAPIT.Blink.Stress_Chlorophyll_A.GWAS.Results.csv', 
            'GAPIT.Blink.Stress_Root_length.GWAS.Results.csv', 
            'GAPIT.Blink.Stress_Root_Na..GWAS.Results.csv', 
            'GAPIT.Blink.Stress_Shoot_K..GWAS.Results.csv', 
            'GAPIT.Blink.Stress_Shoot_length.GWAS.Results.csv', 
            'GAPIT.CMLM.Chlorophyll_A.GWAS.Results.csv', 
            'GAPIT.CMLM.Chlorophyll_B.GWAS.Results.csv', 
            'GAPIT.CMLM.Height_WO_husk.GWAS.Results.csv', 
            'GAPIT.CMLM.Lost_Root_length.GWAS.Results.csv', 
            'GAPIT.CMLM.Lost_Root_thickness.GWAS.Results.csv', 
            'GAPIT.CMLM.Lost_Root_weight.GWAS.Results.csv', 
            'GAPIT.CMLM.Root_K..GWAS.Results.csv', 
            'GAPIT.CMLM.Seed_volume.GWAS.Results.csv', 
            'GAPIT.CMLM.Seed_Weight.GWAS.Results.csv', 
            'GAPIT.CMLM.Shoot_K..GWAS.Results.csv', 
            'GAPIT.CMLM.Shoot_thickness.GWAS.Results.csv', 
            'GAPIT.CMLM.Stress_Root_Na..GWAS.Results.csv', 
            'GAPIT.CMLM.Stress_Shoot_K..GWAS.Results.csv')

makeplots(WF2sug)



makeplotsind <- function(sel) { 
  
  for (file in sel) {
    df <- read_csv(sprintf('../../../GWAS_factory/WF2indNoFDR/%s', file))
    
    name <- paste0(unlist(str_split(file, "\\."))[2:3], collapse = '.')
    
    jpeg(sprintf('%s.jpg', name), width = 1366, height = 788)
    
    manhattan(df, chr="Chromosome", bp="Position", snp="SNP", p="P.value",
              col = c('#9B59B6', '#2E86C1', '#17A589', '#28B463'),
              genomewideline = 6.5,
              annotateTop = True)
    
    dev.off()
  }
}

fv <- str_extract(dir('../../../GWAS_factory/WF2indNoFDR/'), '.*\\Blink.*GWAS.Results.csv')
fv <- fv[!is.na(fv)]

makeplotsind(fv)

indsig <- c('GAPIT.Blink.Seed_volume.GWAS.Results.csv')

makeplotsind(indsig)

indsug <- c('GAPIT.Blink.Lost_Chlorophyll_B.GWAS.Results.csv',
            'GAPIT.Blink.Lost_Root_weight.GWAS.Results.csv',
            'GAPIT.Blink.Lost_Shoot_weight.GWAS.Results.csv',
            'GAPIT.Blink.Seed_volume.GWAS.Results.csv',
            'GAPIT.Blink.Seed_Weight.GWAS.Results.csv',
            'GAPIT.Blink.SES.GWAS.Results.csv',
            'GAPIT.Blink.Shoot_K..GWAS.Results.csv',
            'GAPIT.Blink.Stress_Root_Na..GWAS.Results.csv',
            'GAPIT.Blink.Shoot_thickness.GWAS.Results.csv',
            'GAPIT.Blink.Stress_Chlorophyll_A.GWAS.Results.csv',
            'GAPIT.Blink.Stress_Chlorophyll_B.GWAS.Results.csv',
            'GAPIT.Blink.Stress_Root_K..GWAS.Results.csv',
            'GAPIT.Blink.Stress_Root_thickness.GWAS.Results.csv',
            'GAPIT.Blink.Stress_Shoot_K..GWAS.Results.csv',
            'GAPIT.Blink.Stress_Shoot_length.GWAS.Results.csv',
            'GAPIT.Blink.Stress_Shoot_thickness.GWAS.Results.csv')

makeplotsind(indsug)



makeplotsaus <- function(sel) { 
  
  for (file in sel) {
    df <- read_csv(sprintf('../../../GWAS_factory/WF2ausNoFDR/%s', file))
    
    name <- paste0(unlist(str_split(file, "\\."))[2:3], collapse = '.')
    
    jpeg(sprintf('%s.jpg', name), width = 1366, height = 788)
    
    manhattan(df, chr="Chromosome", bp="Position", snp="SNP", p="P.value",
              col = c('#9B59B6', '#2E86C1', '#17A589', '#28B463'),
              genomewideline = 6.5,
              annotateTop = True)
    
    dev.off()
  }
}

fv <- str_extract(dir('../../../GWAS_factory/WF2ausNoFDR/'), '.*\\Blink.*GWAS.Results.csv')
fv <- fv[!is.na(fv)]

makeplotsaus(fv)


aussig <- c('GAPIT.Blink.Seed_Weight.GWAS.Results.csv')
makeplotsaus(aussig)


setwd('../WF2ausSug')
aussug <- c('GAPIT.Blink.Chlorophyll_A.GWAS.Results.csv',
            'GAPIT.Blink.Chlorophyll_B.GWAS.Results.csv',
            'GAPIT.Blink.Lost_Root_thickness.GWAS.Results.csv',
            'GAPIT.Blink.Lost_Shoot_length.GWAS.Results.csv',
            'GAPIT.Blink.Stress_Chlorophyll_A.GWAS.Results.csv',
            'GAPIT.Blink.Stress_Root_K..GWAS.Results.csv',
            'GAPIT.Blink.Stress_Shoot_K..GWAS.Results.csv',
            'GAPIT.Blink.Width_WO_husk.GWAS.Results.csv',
            'GAPIT.CMLM.Width_WO_husk.GWAS.Results.csv')

makeplotsaus(aussug)


#make qq
library(qqman)
library(tidyverse)

makeplots <- function(sel) { 
  
  for (file in sel) {
    df <- read_csv(sprintf('../../../GWAS_factory/WF2noFDR/%s', file))
    
    name <- paste0(unlist(str_split(file, "\\."))[2:3], collapse = '.')
    
    jpeg(sprintf('%s.jpg', name), width = 788, height = 788)
    
    qq(df$P.value)
    
    dev.off()
  }
}

fv <- str_extract(dir('../../../GWAS_factory/WF2noFDR/'), '.*\\Blink.*GWAS.Results.csv')
fv <- fv[!is.na(fv)]

makeplots(fv)


makeplots <- function(sel) { 
  
  for (file in sel) {
    df <- read_csv(sprintf('../../../GWAS_factory/WF2ausNoFDR/%s', file))
    
    name <- paste0(unlist(str_split(file, "\\."))[2:3], collapse = '.')
    
    jpeg(sprintf('%s.jpg', name), width = 788, height = 788)
    
    qq(df$P.value)
    
    dev.off()
  }
}

fv <- str_extract(dir('../../../GWAS_factory/WF2ausNoFDR/'), '.*\\Blink.*GWAS.Results.csv')
fv <- fv[!is.na(fv)]

makeplots(fv)



makeplots <- function(sel) { 
  
  for (file in sel) {
    df <- read_csv(sprintf('../../../GWAS_factory/WF2indNoFDR/%s', file))
    
    name <- paste0(unlist(str_split(file, "\\."))[2:3], collapse = '.')
    
    jpeg(sprintf('%s.jpg', name), width = 788, height = 788)
    
    qq(df$P.value)
    
    dev.off()
  }
}

fv <- str_extract(dir('../../../GWAS_factory/WF2indNoFDR/'), '.*\\Blink.*GWAS.Results.csv')
fv <- fv[!is.na(fv)]

makeplots(fv)
