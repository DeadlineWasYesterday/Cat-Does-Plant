
ou <- read.table('../Data/BLASTout/saltyexnpbrout.csv',
                sep = '\t',
               col.names = c('q','acc','pi','len','mm','gap','qs','qe','ss','se','e','sc'))

ou <- merge(ou, select(data.frame(table(ou$q)), q = Var1, Freq),
            by = 'q', all = TRUE)
ou <- ou[order(ou$pi, decreasing = TRUE),]

