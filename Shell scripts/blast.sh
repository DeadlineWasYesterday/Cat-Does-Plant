#!/bin/bash

makeblastdb -in IRGSP-1.0_genome.fasta -dbtype nucl

blastn -query IRGSP-1.0_gene_2020-09-09.fasta -db IRGSP-1.0_genome.fasta -qcov_hsp_perc 100 -outfmt 6 -out BoutAllgenes.tsv

-query RAPDB\ all\ exons.fasta -db IRGSP-1.0_genome.fasta -qcov_hsp_perc 100 -outfmt 6 -out BoutAllexons.tsv