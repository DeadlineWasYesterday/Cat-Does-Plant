#!/bin/bash

for f in *
do
zcat $f | bgzip -c --threads 30 > ../../Cat-chan/176vcfbz/$f && tabix ../../Cat-chan/176vcfbz/$f
done
