# DMS2structure

These are scripts used for the analysis of deep mutational scanning data in Schmiedel & Lehner, bioRxiv 2018 (https://www.biorxiv.org/content/early/2018/04/20/303875)  

The DATASET_pipeline.R scripts do the complete analysis for one dataset, the necessary data are already deposited in the [dataset_folder]/dataset/ folder, except for the WW domain, for which the sequencing data has to be downloaded and processed separately (as described in the pipeline script).  


## required software
To run these scripts, you will need   
R (version 3.4)  
DeepContact (can be cloned from https://github.com/largelymfs/deepcontact )  
XPLOR-NIH (can be downloaded from https://nmr.cit.nih.gov/xplor-nih/ )  
TMscore (can be downloaded from https://zhanglab.ccmb.med.umich.edu/TM-score/ )  

## required R packages
data.table  
ggplot2  
cowplot  
GGally  
mgcv  
caTools  
parallel  
stringr  
gdata  
corpcor  
Rpdb  
pdist  
metap  
RColorBrewer  
ssh.utils  
seqinr  
optparse  
pheatmap