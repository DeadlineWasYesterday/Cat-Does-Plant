
mtbl <- read_csv('../Potential Salt Genes/gtbl 19cols.csv', col_types = cols(.default = col_character()))
fnames <- paste(mtbl$chr, mtbl$withpro, mtbl$end, mtbl$wGName, sep = "-")

fnames <- lapply(fnames, function(name) { 
  name <- gsub('\\|', '-', name)
  name <- gsub('\\*', '(star)', name)
  name <- gsub('\\/', '(slash)', name)
  name <- gsub('\\:', '(colon)', name)
  #name <- gsub('\\', '(slash2)', name)
} )

genes <- read.fasta('../Pie/884 salty genes.fasta', seqtype = "DNA", as.string = FALSE, forceDNAtolower = FALSE, whole.header = FALSE)

for (name in fnames) {
  
  snp <- read_csv(sprintf('../Potential Salt Genes/snpia/%s.csv', name), col_types = cols(.default = col_character()))
  snp <- snp[1,]
  
  lcs <- (unlist(lapply(colnames(snp), strtoi)))
  lcs <- lcs[!is.na(lcs)]
  lcs <- sapply(lcs, paste)
  
  snp <- snp[lcs]
  
  wpstart <- strtoi(strsplit(name, '-')[[1]][2])
  gname <- strsplit(name, '^[A-Z0-9]*-*[0-9]*-*[0-9]*-')[[1]][2]
  
  gname <- gsub("(star)", "*", gname, fixed = TRUE)
  gname <- gsub("(slash)", "/", gname, fixed = TRUE)
  gname <- gsub("(colon)", ":", gname, fixed = TRUE)
  
  start <- mtbl$start[mtbl$wGName == gname]
  lcid <- mtbl$tr1[mtbl$wGName == gname]
  lcid <- strsplit(lcid, '\\.')[[1]][1]
  
  seq <- genes[lcid]
  seq <- seq[[1]]
  
  posvecseq <- unlist(lapply(colnames(snp), strtoi)) - strtoi(start)
  posvecseq <- posvecseq[which(posvecseq > 0)] + 1
  
  posveccol <- colnames(snp)[which(colnames(snp) > start)]
  
  n <- 0
  if (all(seq[posvecseq] == snp[posveccol])) {
    mtbl[which(mtbl$wGName == gname), 'QC'] <- paste(0)
  }
  else {
  for (i in 1:500) {
    
    upcon <- posvecseq + i
    upcon <- upcon[which(upcon > 0)]
    downcon <- posvecseq + i
    downcon <- downcon[which(downcon > 0)]
    
    if (all(seq[upcon] == snp[posveccol])) {
      n <- i
      break 
      #subtract n from exon position 
      }
      else if (all(seq[downcon] == snp[posveccol])) {
        n <- -i
        break
      #add n to exon position
      }
    }
  }
    mtbl[which(mtbl$wGName == gname), 'QC'] <- paste(n)
}
