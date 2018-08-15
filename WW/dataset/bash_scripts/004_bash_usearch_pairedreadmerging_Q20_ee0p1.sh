#!/bin/bash
#$ -q long-sl7
#$ -N usearch_prm
#$ -pe smp 4
#$ -l virtual_free=16G,h_rt=42:00:00
#$ -o DMS2struct/datasets/WW_Araya2012
#$ -e DMS2struct/datasets/WW_Araya2012

# command to run usearch paired read merger
/software/bl/el7.2/usearch-10.0.240/usearch10.0.240_i86linux32 -fastq_mergepairs /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ/SRR569005_R1.fastq -fastqout /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Input1_Q20ee0p1.fastq -report /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Input1_Q20ee0p1.report -fastq_minqual 20 -fastq_merge_maxee 0.01 

/software/bl/el7.2/usearch-10.0.240/usearch10.0.240_i86linux32 -fastq_mergepairs /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ/SRR569006_R1.fastq -fastqout /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Input2_Q20ee0p1.fastq -report /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Input2_Q20ee0p1.report -fastq_minqual 20 -fastq_merge_maxee 0.01 

/software/bl/el7.2/usearch-10.0.240/usearch10.0.240_i86linux32 -fastq_mergepairs /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ/SRR569007_R1.fastq -fastqout /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round1-1_Q20ee0p1.fastq -report /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round1-1_Q20ee0p1.report -fastq_minqual 20 -fastq_merge_maxee 0.01 

/software/bl/el7.2/usearch-10.0.240/usearch10.0.240_i86linux32 -fastq_mergepairs /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ/SRR569008_R1.fastq -fastqout /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round1-2_Q20ee0p1.fastq -report /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round1-2_Q20ee0p1.report -fastq_minqual 20 -fastq_merge_maxee 0.01 

/software/bl/el7.2/usearch-10.0.240/usearch10.0.240_i86linux32 -fastq_mergepairs /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ/SRR569009_R1.fastq -fastqout /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round2-1_Q20ee0p1.fastq -report /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round2-1_Q20ee0p1.report -fastq_minqual 20 -fastq_merge_maxee 0.01 

/software/bl/el7.2/usearch-10.0.240/usearch10.0.240_i86linux32 -fastq_mergepairs /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ/SRR569010_R1.fastq -fastqout /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round2-2_Q20ee0p1.fastq -report /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round2-2_Q20ee0p1.report -fastq_minqual 20 -fastq_merge_maxee 0.01 

/software/bl/el7.2/usearch-10.0.240/usearch10.0.240_i86linux32 -fastq_mergepairs /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ/SRR569011_R1.fastq -fastqout /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round3-1_Q20ee0p1.fastq -report /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round3-1_Q20ee0p1.report -fastq_minqual 20 -fastq_merge_maxee 0.01 

/software/bl/el7.2/usearch-10.0.240/usearch10.0.240_i86linux32 -fastq_mergepairs /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ/SRR569012_R1.fastq -fastqout /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round3-2_Q20ee0p1.fastq -report /users/blehner/jschmiedel/DMS2struct/datasets/WW_Araya2012/FastQ_merged/WW1_Round3-2_Q20ee0p1.report -fastq_minqual 20 -fastq_merge_maxee 0.01 
