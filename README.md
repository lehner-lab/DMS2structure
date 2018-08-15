# DMS2structure

These are scripts used for the analysis of deep mutational scanning data in Schmiedel & Lehner, bioRxiv 2018 (https://www.biorxiv.org/content/early/2018/04/20/303875)\n

The DATASET_pipeline.R scripts do the complete analysis for one dataset, the necessary data are already deposited in the [dataset_folder]/dataset/ folder, except for the WW domain, for which the sequencing data has to be downloaded and processed separately (as described in the pipeline script).\n


## required software
To run these scripts, you will need\n 
R (version 3.4)\n
DeepContact (can be cloned from https://github.com/largelymfs/deepcontact )\n
XPLOR-NIH (can be downloaded from https://nmr.cit.nih.gov/xplor-nih/ )\n
TMscore (can be downloaded from https://zhanglab.ccmb.med.umich.edu/TM-score/ )\n

## required R packages
data.table\n
ggplot2\n
cowplot\n
GGally\n
mgcv\n
caTools\n
parallel\n
stringr\n
gdata\n
corpcor\n
Rpdb\n
pdist\n
metap\n
RColorBrewer\n
ssh.utils\n
seqinr\n
optparse\n
pheatmap\n
