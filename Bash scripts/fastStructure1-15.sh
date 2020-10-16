#!/bin/bash


#for whole population
python structure.py -K 1 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 2 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 3 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 4 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 5 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 6 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 7 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 8 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 9 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 10 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 11 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 12 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 13 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 14 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs
python structure.py -K 15 --input=../../Data/WF2plinkbed --output=../fsout/WF2fs



python chooseK.py --input=../fsout/WF2fs

#Model complexity that maximizes marginal likelihood = 3
#Model components used to explain structure in data = 4


python distruct.py -K 3 --input=../fsout/WF2fs --output=../fsout/fig --title='All_accessions'