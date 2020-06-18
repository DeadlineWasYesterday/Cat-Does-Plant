
mtbl <- read_csv('../Potential Salt Genes/gtbl 19cols.csv', col_types = cols(.default = col_character()))
fnames <- paste(mtbl$chr, mtbl$withpro, mtbl$end, mtbl$wGName, sep = "-")

mtbl <- mtbl[-633,] #duplicate

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
  end <- strtoi(strsplit(name, '-')[[1]][3])
  start <- strtoi(wpstart) + 2000
  
  tr1 <- mtbl$tr1[mtbl$end == end]
  lcid <- strsplit(tr1, '\\.')[[1]][1]
  
  seq <- genes[lcid]
  if (is.null(seq)) {
    mtbl[which(mtbl$tr1 == tr1), 'QC'] <- 'miss'
    next
  }
  if (name == "CHR7-20293862-20297855-OsGPT2") {
    mtbl[which(mtbl$tr1 == tr1), 'QC'] <- 'miss'
    next
  }
  seq <- seq[[1]]
  
  posvecseq <- unlist(lapply(colnames(snp), strtoi)) - strtoi(start)
  posvecseq <- posvecseq + 1
  posvecseq <- posvecseq[which(posvecseq > 0)]
  
  posveccol <- colnames(snp)[which(colnames(snp) > start)]
  
  n <- 'def'
  if (length(posvecseq) != length(posveccol) && 
      length(posveccol) >= 10) { 
     posveccol <- posveccol[1:10]
     posvecseq <- posvecseq[1:10]
  }
  if (any(is.na(seq[posvecseq])) || any(is.na(snp[posveccol]))) {
    posveccol <- posveccol[1:10]
    posvecseq <- posvecseq[1:10]
if (any(is.na(seq[posvecseq])) || any(is.na(snp[posveccol]))){n<-'fail'}
  }
  
  
  if (all(seq[posvecseq] == snp[posveccol])) {
    mtbl[which(mtbl$tr1 == tr1), 'QC'] <- 'pass'
    n <- 'pass'
  }
  
  else {
  for (i in 1:500) {
    upcon <- posvecseq + i
    upcon <- upcon[which(upcon > 0)]
    downcon <- posvecseq - i
    downcon <- downcon[which(downcon > 0)]
    
    if (any(is.na(seq[upcon])) || any(is.na(seq[downcon])) 
        || length(upcon) == 0 || length(downcon) == 0) {
      posveccol <- posveccol[1:10]
      posvecseq <- posvecseq[1:10]
      upcon <- upcon[1:10]
      downcon <- downcon[1:10]
      if (any(is.na(seq[upcon])) || any(is.na(seq[downcon])) 
          || length(upcon) == 0 || length(downcon) == 0) {
      n <- 'fail'
      break }
    }
    
    if (all(seq[upcon] == snp[posveccol])) {
      n <- i
      break 
      #subtract n from exon position 
    }
    
    if (all(seq[downcon] == snp[posveccol])) {
      n <- -i
      break
      #add n to exon position
      }
  }
    mtbl[which(mtbl$tr1 == tr1), 'QC'] <- paste(n)
  }

  } 
