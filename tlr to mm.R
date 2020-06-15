library(tidyverse)

for (file in dir('../Potential Salt Genes/snpout/')) {
  snpsc <- read_csv(sprintf('../Potential Salt Genes/snpout/%s', file))
  
  png(filename=sprintf('../Potential Salt Genes/plots/sesvsmm%s.png', file))
  plot(snpsc$SES, snpsc$MISMATCH, xlab = 'SES', 
       ylab = 'Mismatches', main = file)
  
  abline(lm(snpsc$MISMATCH ~ snpsc$SES, data = snpsc), col = "blue")
  
  
  dev.off()
  
}

a = c(8,7,6,4,2)
b = c(2,3,12,17,25)
tab <- data.frame(a, b)

plot(a,b)
abline(lm(b ~ a, data = tab), col = "blue")

slope <- lm(b ~ a, data = tab)$coefficients[2]