library(devtools)
install_github("isglobal-brge/SNPassoc")
library(SNPassoc)

# Next steps may be very time consuming. So they are not executed
data(HapMap)
myDat<-setupSNP(HapMap, colSNPs=3:1000, sort = TRUE,
 info=HapMap.SNPs.pos, sep="")

resHapMap<-scanWGassociation(group~1, data=myDat, model="log", nperm = 200)

permTest(resHapMap, 'minimum')
