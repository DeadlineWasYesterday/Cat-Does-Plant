#!/bin/bash

./PopLDdecay -InVCF ../Data/WF2indR7.vcf -MAF 0.2 -OutStat WF2indR7
./Plot_OnePop.pl -InFile WF2indR7.stat.gz -output WF2indR7


./PopLDdecay -InVCF ../Data/WF2ausR7.vcf -MAF 0.2 -OutStat WF2ausR7
./Plot_OnePop.pl -InFile WF2ausR7.stat.gz -output WF2ausR7

./PopLDdecay -InVCF ../Data/WF1M20.vcf -MAF 0.25 -OutStat WF1M20
./Plot_OnePop.pl -InFile WF1M20.stat.gz -output WF1M20


#these are good notes below, but we did not use them.
#from https://www.biostars.org/p/300381/#300423

## plink --bfile "snp" --r2 --out "snp"

## mapthin -b 20 snp.bim snp-thin.bim

## plink --bfile "snp-thin" --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out "snp-thin"

## cat snp-thin.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > snp-thin.ld.summary
