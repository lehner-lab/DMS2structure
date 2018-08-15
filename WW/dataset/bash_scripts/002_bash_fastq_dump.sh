#!/bin/bash
# -q long-sl7
# -N sra.download
# -l virtual_free=10G,h_rt=42:00:00
# -o DMS2struct/datasets/WW_Araya2012
# -e DMS2struct/datasets/WW_Araya2012

# go to directory where all the sra files are located
cd DMS2struct/datasets/WW_Araya2012
mkdir FastQ/

# run the code on all the .sra files I downloaded just now.
for filename in ./*.sra

do
    /users/blehner/pbaeza/bin/fastq-dump --outdir cd /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ --split-files ${filename}

done

