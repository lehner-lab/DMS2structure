#!/bin/bash
# -q long-sl7
# -N fastqc
# -l virtual_free=10G,h_rt=42:00:00
# -o DMS2struct/datasets/WW_Araya2012
# -e DMS2struct/datasets/WW_Araya2012

# command to run FastQC
# - o argument is output folder which previously exists
mkdir DMS2struct/datasets/WW_Araya2012/FastQC/
# last argument is the file to be processed by fastqc
/software/bl/el6.5/fastqc/FastQC/fastqc -o /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQC/ /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ/*


