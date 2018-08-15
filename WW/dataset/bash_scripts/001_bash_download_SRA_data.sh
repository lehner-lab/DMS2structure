#!/bin/bash
# -q short-sl7,long-sl7
# -N sra.download
# -l virtual_free=1G,h_rt=42:00:00
# -o DMS2struct/datasets/WW_Araya2012
# -e DMS2struct/datasets/WW_Araya2012

### download sequencing data from Araya 2012

mkdir DMS2struct/datasets/WW_Araya2012/
cd DMS2struct/datasets/WW_Araya2012/

# download all SRR files
wget -r -nH -nd -np -R index.html* ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP015/SRP015751/

