#!/bin/bash
#$ -q long-sl7
#$ -N usearch_prm
#$ -pe smp 4
#$ -l virtual_free=16G
#$ -M joern.schmiedel@crg.eu
#$ -m ae
#$ -o DMS2struct/datasets/WW_Araya2012
#$ -e DMS2struct/datasets/WW_Araya2012

# command to run usearch fastx_unique code to collapse unique variants

for filename in /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/*Q20ee0p1.fastq

do
	/software/bl/el7.2/usearch-10.0.240/usearch10.0.240_i86linux32 -fastx_uniques ${filename} -fastaout ${filename}_unique -sizeout -relabel Uniq

done
