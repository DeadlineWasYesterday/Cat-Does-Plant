#!/bin/bash


vcftools --vcf WF2.vcf --maf 0.05 --plink --out WF2


#After filtering, kept 176 out of 176 Individuals
#Writing PLINK PED and MAP files ... 
#Done.
#After filtering, kept 4364394 out of a possible 9260197 Sites
#Run Time = 405.00 seconds

plink --map WF2.map --ped WF2.ped --maf 0.05 --indep-pairwise 50 50 .2 --out WF2

#Use WF2.pruned.in to generate WF2LE.vcf

vcftools --vcf WF2LE.vcf --maf 0.05 --plink --out WF2LE

plink --file WF2LE --make-bed --out WF2LE





#Pruning WF5R8
vcftools --vcf WF5R8.vcf --maf 0.05 --plink --out WF5R8
plink --map WF5R8.map --ped WF5R8.ped --maf 0.05 --indep-pairwise 50 50 .2 --out WF5R8
#Use WF5R8.pruned.in to generate WF5R8LE.vcf

