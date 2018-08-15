####################################
##### WW domain from Araya2012 #####
####################################


#this is the pipeline used to analyse the WW domain data from Araya et al. 2012

#first, set the working directory to the DMS2structure folder
setwd("/where/is/DMS2structure/")


#source scripts
filelist = list.files('scripts/')
sapply(paste0('scripts/',filelist),source,.GlobalEnv)

#create the necessary subfolder structure for all results and processed data
dataset_dir = "WW/"
create_directory_structure(dataset_dir)
#then save this script in the dataset_dir
#and paste all necessary source data into dataset_dir/dataset

#load required packages
require(data.table)
require(ggplot2)
require(cowplot)
require(GGally)
require(seqinr)
theme_set(theme_minimal())



###################################################
##### download and preprocess sequencing data #####
###################################################

#to download data from the Sequencing Read Archive, execute the bash scripts (via qsub) in the dataset/bash_script/ folder
#you will need wget, fastq-dumb, fastqc (optional) and usearch programs
#(and adjust the folder/file pointers in the bash scripts)
#copy the final output files after merging and unique variant counting to the dataset_dir/dataset/sequencing_data/ directory

#### load usearch mapped sequencing reads [takes long, skip to section below to load data]
filenames = dir(path=paste0(dataset_dir,"/dataset/sequencing_data/"))
filenames = filenames[grep(pattern = "*.fastq_unique",x = filenames)]

file_id = sapply(X=filenames,FUN = function(X) {strsplit(X,"_")[[1]][2]})

wt_NTseq = tolower("GACGTTCCACTGCCGGCTGGTTGGGAAATGGCTAAAACTAGTTCTGGTCAGCGTTACTTCCTGAACCACATCGACCAGACCACCACGTGGCAGGACCCGCGT")
wt_NTseq_split = strsplit(wt_NTseq,"")[[1]]
wt_AAseq = paste0(seqinr::translate(strsplit(wt_NTseq,split="")[[1]]),collapse="")
write(paste0(">WW domain aaseq\n",wt_AAseq),file = paste0(dataset_dir,"dataset/WW_sequence.fasta"))
wt_AAseq_split = strsplit(wt_AAseq,"")[[1]]

for (i in 1:length(filenames)) {
  print(i)
  output = seqinr::read.fasta(file = paste0(dataset_dir,"/dataset/sequencing_data/",filenames[i]),as.string = T)
  #read nucleotide sequence
  nt_seq = sapply(X=1:length(output), FUN = function(X){as.character(output[[X]])})
  #read number of appearances in the sequencing file
  count = sapply(X=1:length(output), FUN = function(X){as.numeric(strsplit(attr(output[[X]],"name"),split="[;=]")[[1]][3])})

  fq = data.table(nt_seq,count,length=nchar(nt_seq))
  fq = fq[length == 102] #only keep sequences with correct length
  fq[,aa_seq := paste0(seqinr::translate(strsplit(nt_seq,split="")[[1]]),collapse=""),nt_seq] #translate to AA seq
  fq[,Nmut_aa := sum(strsplit(aa_seq,"")[[1]] != wt_AAseq_split),aa_seq] #count AA mutations
  fq = fq[Nmut_aa <= 2] #only keep variants with at max 2 AA mutations
  names(fq)[2] = paste0("count_", file_id[i])
  
  #merge data from different replicates/timepoints
  if (i==1) {
    variant_data = fq
  } else {
    variant_data = merge(variant_data,fq[,.(nt_seq,.SD),,.SDcols=names(fq)[grep(names(fq),pattern="count")]],by="nt_seq",all = T)
    names(variant_data)[grep(names(variant_data),pattern=".SD")] = paste0("count_", file_id[i])
  }
}
names(variant_data) = c("nt_seq","count_r1_t0","length","aa_seq","Nmut_aa","count_r2_t0","count_r1_t1","count_r2_t1","count_r1_t2",
                        "count_r2_t2","count_r1_t3","count_r2_t3")

#fill in aa_seq & Nmut_aa for variants only contained in one of the later files
variant_data[is.na(aa_seq),aa_seq := paste0(seqinr::translate(strsplit(nt_seq,split="")[[1]]),collapse=""),nt_seq]
variant_data[is.na(Nmut_aa),Nmut_aa := sum(strsplit(aa_seq,"")[[1]] != wt_AAseq_split),aa_seq]

#also count nucleotide mutations for all variants
variant_data[,Nmut_nt := sum(strsplit(nt_seq,"")[[1]] != wt_NTseq_split),nt_seq]

# save the processed data
write.table(x = variant_data,file = paste0(dataset_dir,"processed_data/mapped_variant_data_Q20.txt"),
            row.names = F,col.names = T,quote = F)




#############################################################################################
##### preprocess data (calculate fitness scores and errors, set quality thresholds etc) #####
#############################################################################################

###### load preprocessed data (if repeating the script)
variant_data = fread(paste0(dataset_dir,"processed_data/mapped_variant_data_Q20.txt"))


## statistics on variant mutations and total read counts
variant_data[,.(.N,count_r1 = sum(count_r1_t0,na.rm=T),
                count_r2 = sum(count_r2_t0,na.rm=T)),.(Nmut_aa,Nmut_nt)][order(Nmut_aa,Nmut_nt)]

###### merge synonymous variants
## use double mutants + 1nt away
variant_data0 = unique(variant_data[!is.na(aa_seq) & Nmut_aa >= (Nmut_nt-1),
                             .(Nmut_aa,count_r1_t0 = sum(count_r1_t0,na.rm=T),
                               count_r1_t1 = sum(count_r1_t1,na.rm=T),
                               count_r1_t2 = sum(count_r1_t2,na.rm=T),
                               count_r1_t3 = sum(count_r1_t3,na.rm=T),
                               count_r2_t0 = sum(count_r2_t0,na.rm=T),
                               count_r2_t1 = sum(count_r2_t1,na.rm=T),
                               count_r2_t2 = sum(count_r2_t2,na.rm=T),
                               count_r2_t3 = sum(count_r2_t3,na.rm=T)),aa_seq])

## only retain variants properly observed in one of the two replicatds
variant_data1 = variant_data0[(count_r1_t0 > 10 & count_r1_t1 > 0 ) | (count_r2_t0 > 10  & count_r2_t1 > 0)]

## mark the wild-type
variant_data1[aa_seq == wt_AAseq,WT := TRUE]
## mark variants with a STOP codon (marked by *)
variant_data1[,STOP := ifelse(length(grep(aa_seq,pattern="\\*"))==1,TRUE,FALSE),aa_seq]

###### calculate fitness and error of estimates over time-points
dataset = copy(variant_data1)
TP = 4
Nrep = 2
for (r in 1:Nrep) {
  # define X to calculate slope
  # X = matrix(rep(0:(TP-1),each=nrow(dataset)),nrow=nrow(dataset),ncol=TP)
  #spacing is not equidistant
  X = matrix(rep(cumsum(c(0,0.6,1.17,1.22)),each=nrow(dataset)),nrow=nrow(dataset),ncol=TP)
  
  for (t in 0:(TP-1)) {
    #time-point/replicate scores
    dataset[,paste0("score_r",r,"_t",t) := log(get(paste0("count_r",r,"_t",t))/dataset[WT == T,get(paste0("count_r",r,"_t",t))])]
    #time-point/replicate error
    dataset[,paste0("sigma_r",r,"_t",t) := sqrt(1/get(paste0("count_r",r,"_t",t)) + 1/dataset[WT == T,get(paste0("count_r",r,"_t",t))])]
  }
  #weighted slopes (Barlow1989 section 6.2 p101f)
  score_rx_tx = dataset[,.SD,.SDcols = grep(paste0("score_r",r,"_t[0-9]"),colnames(dataset))]
  sigma_rx_tx = dataset[,.SD,.SDcols = grep(paste0("sigma_r",r,"_t[0-9]"),colnames(dataset))]
  
  X[is.na(score_rx_tx) | score_rx_tx == -Inf] = NA
  sigma_rx_tx[is.na(score_rx_tx) | score_rx_tx == -Inf] = NA
  score_rx_tx[score_rx_tx == -Inf] = NA
  #weighted slopes
  dataset[,paste0("score_r",r) := (rowSums(score_rx_tx*X/sigma_rx_tx^2,na.rm=T)/rowSums(1/sigma_rx_tx^2,na.rm=T) - 
                                     rowSums(X/sigma_rx_tx^2,na.rm=T)/rowSums(1/sigma_rx_tx^2,na.rm=T) * 
                                     rowSums(score_rx_tx/sigma_rx_tx^2,na.rm=T)/rowSums(1/sigma_rx_tx^2,na.rm=T)) /
            (rowSums(X^2/sigma_rx_tx^2,na.rm=T)/rowSums(1/sigma_rx_tx^2,na.rm=T) - 
               (rowSums(X/sigma_rx_tx^2,na.rm=T)/rowSums(1/sigma_rx_tx^2,na.rm=T))^2)]
  #weighted errors (Barlow1989 section 6.2 p101f)
  dataset[,paste0("sigma_r",r) := sqrt( 1/rowSums(1/sigma_rx_tx^2,na.rm=T) / 
                                          (rowSums(X^2/sigma_rx_tx^2,na.rm=T)/rowSums(1/sigma_rx_tx^2,na.rm=T) -
                                             (rowSums(X/sigma_rx_tx^2,na.rm=T)/rowSums(1/sigma_rx_tx^2,na.rm=T))^2))]
}

# plot replicate scores against each other
ggplot(dataset,aes(score_r1,score_r2)) +
  geom_hex()

#get rid of few variants with way too low scores in both replicates
dataset = dataset[score_r1 > -1.5 & score_r2 > -1.5] 


###### merge replicates
score_rx = dataset[,.SD,.SDcols = grep(paste0("score_r[0-9]$"),colnames(dataset))]
sigma_rx = dataset[,.SD,.SDcols = grep(paste0("sigma_r[0-9]$"),colnames(dataset))]
dataset[,fitness := rowSums(score_rx/(sigma_rx^2))/rowSums(1/(sigma_rx^2))]
dataset[,sigma := sqrt(1/rowSums(1/(sigma_rx^2)))]

# plot fitness distribution
ggplot(dataset,
       aes(fitness,color=factor(Nmut_aa))) +
  geom_density() +
  theme_minimal() +
  labs(x = "log(fitness)",color = "#mut")
ggsave(paste0(dataset_dir,"results/preprocessing/fitness_distribution.pdf"),width=4,height = 3)

# plot fitnes versus error
ggplot(dataset[Nmut_aa > 0],
       aes(fitness,sigma)) +
  geom_hex() +
  scale_fill_distiller(direction=1,trans="log10") +
  scale_y_log10() +
  facet_wrap( ~ Nmut_aa) +
  theme_minimal() +
  labs(x = "log(fitness)",y = "error of fitness estimate")
ggsave(paste0(dataset_dir,"results/preprocessing/fitness_error_distribution.pdf"),width=8,height=4)



###### build wildtype, singles and doubles data.tables
wildtype = dataset[Nmut_aa==0,.(fitness,sigma,STOP)]

singles = dataset[Nmut_aa==1,.(mean_count_t0 = mean(c(count_r1_t0,count_r2_t0)),
                               fitness,sigma,STOP),aa_seq]
singles[,Pos := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split),aa_seq]
singles[,Mut := strsplit(aa_seq,"")[[1]][Pos],aa_seq]
singles[,WT_AA := wt_AAseq_split[Pos],aa_seq]
singles[,aa_seq := NULL]
setcolorder(singles,names(singles)[c(5:7,1:4)])

doubles = dataset[Nmut_aa==2,.(mean_count_t0 = mean(c(count_r1_t0,count_r2_t0)),fitness,sigma,STOP),aa_seq]
doubles[,Pos1 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[1],aa_seq]
doubles[,Pos2 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[2],aa_seq]
doubles[,Mut1 := strsplit(aa_seq,"")[[1]][Pos1],aa_seq]
doubles[,Mut2 := strsplit(aa_seq,"")[[1]][Pos2],aa_seq]
doubles[,WT_AA1 := wt_AAseq_split[Pos1],aa_seq]
doubles[,WT_AA2 := wt_AAseq_split[Pos2],aa_seq]
doubles[,aa_seq := NULL]

###### pull single mutant fitness and error into double mutant table
doubles[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],by=.(Pos1,Mut1)]
doubles[,sigma1 := singles[Pos == Pos1 & Mut == Mut1,sigma],by=.(Pos1,Mut1)]
doubles[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],by=.(Pos2,Mut2)]
doubles[,sigma2 := singles[Pos == Pos2 & Mut == Mut2,sigma],by=.(Pos2,Mut2)]


###### mark variants with nonsensical fitness values / STOP variants
wildtype[,is.fitness := TRUE]
singles[,is.fitness := !STOP]
doubles[,is.fitness := !STOP]

###### estimate lower measurement limit of fitness assay 
# from STOP variants
ggplot(rbind(singles[,.(fitness,STOP,type="singles")],
             doubles[,.(fitness,STOP,type="doubles")]),
       aes(fitness,color=STOP,linetype = type)) + 
  geom_density() +
  geom_vline(xintercept = singles[STOP==T,weighted.mean(fitness,sigma^-2)],linetype=2) +
  geom_vline(xintercept = doubles[STOP==T,weighted.mean(fitness,sigma^-2)],linetype=1)
ggsave(paste0(dataset_dir,"results/preprocessing/lower_fitnes_bound.pdf"),width=5,height=4)
#same for singles and doubles
lower_bound_F = singles[STOP==T,weighted.mean(fitness,sigma^-2)]
print(lower_bound_F)
doubles[STOP==T,weighted.mean(fitness,sigma^-2)]

###### define which variants have enough read coverage
# plot read distribution
ggplot(dataset[Nmut_aa>0,.(mean_count_t0 = (count_r1_t0+count_r2_t0)/2,fitness,Nmut_aa)],
       aes(mean_count_t0,fitness,color=factor(Nmut_aa))) +
  geom_hex() +
  scale_x_log10() +
  theme_minimal() +
  labs(x = "read counts in input library",
       y = "read counts in output library")
ggsave(paste0(dataset_dir,"results/preprocessing/readcounts_input_versus_fitness.pdf"),width=6,height=6)
# > there doesn't seem to be any low input coverage limitations, discard few doubles with <= 10 average input reads
wildtype[,is.reads0 := TRUE]
singles[,is.reads0 := TRUE]
doubles[mean_count_t0 > 10,is.reads0 := TRUE]


###### reorder doubles table
doubles = doubles[,.(Pos1,WT_AA1,Mut1,Pos2,WT_AA2,Mut2,
                    STOP,mean_count_t0,is.fitness,is.reads0,
                    fitness1,sigma1,fitness2,sigma2,
                    fitness,sigma)]

###### write data to files
write.table(x = wildtype, file = paste0(dataset_dir,"processed_data/DMS_wildtype.txt"),
            quote = F,row.names = F, col.names = T)
write.table(x = singles, file = paste0(dataset_dir,"processed_data/DMS_singles.txt"),
            quote = F,row.names = F, col.names = T)
write.table(x = doubles, file = paste0(dataset_dir,"processed_data/DMS_doubles_preE.txt"),
            quote = F,row.names = F, col.names = T)





##################################################################
### calculate epistasis null model and call pos./neg.epistasis ###
##################################################################

prefix = "WW_"

#call epistatic interactions
# doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles_preE.txt"))
doubles = call_epistasis_binary(doubles,
                                lower_bound_F,
                                dataset_dir = dataset_dir,
                                prefix = prefix)


#############################################
### calculate pairwise interaction scores ###
#############################################

# doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles.txt"))
PWI = calculate_pairwise_interaction_scores(doubles,
                                            dataset_dir = dataset_dir,
                                            detailed = T)


#########################################
### deep contact transform PWI scores ###
#########################################

##### full range WW domain
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))
PWI_transformed = deepcontact_transform_basic2d(PWI[,.(Pos1,Pos2,
                                                       epistasis_score,association_score,combined_score)],
                                                dataset_dir = dataset_dir,
                                                deepcontact_dir = "where/is/deepcontact/",
                                                prefix = prefix)

### negative controls for DeepContact
# 3x permutate combined_scores, while keeping matrix symmetry
set.seed(1603)
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))[Pos1<Pos2,.(Pos1,Pos2,control1 = sample(combined_score),
                                                                          control2 = sample(combined_score),
                                                                          control3 = sample(combined_score))]
PWI2 = switch_double_DT(PWI,list(c("Pos1","Pos2")),c("control1","control2","control3"))
PWI2 = rbind(PWI2,data.table(Pos1=1:34,Pos2=1:34,control1=NA,control2=NA,control3=NA))
write.table(x = PWI2, file = paste0(dataset_dir,"processed_data/DMS_PWI_DC_control.txt"),
            quote = F,row.names = F, col.names = T)
PWI_transformed = deepcontact_transform_basic2d(PWI2,
                                                dataset_dir = dataset_dir,
                                                output_filename = "DMS_PWI_DC_control_deepcontact.txt",
                                                deepcontact_dir = "where/is/deepcontact/",
                                                prefix = prefix)

##### do a DeepContact transform using only positions 1-29 to avoid artifacts
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))
PWI = PWI[Pos1 < 30 & Pos2 < 30]
write.table(x = PWI, file = paste0(dataset_dir,"/processed_data/DMS_PWI_1to29.txt"),
            quote = F,row.names = F, col.names = T)
PWI_transformed = deepcontact_transform_basic2d(PWI[,.(Pos1,Pos2,
                                                       epistasis_score,association_score,combined_score)],
                                                dataset_dir = dataset_dir,
                                                output_filename = "DMS_PWI_1to29_deepcontact.txt",
                                                deepcontact_dir = "where/is/deepcontact/",
                                                prefix = paste0(prefix,"1to29_"))


### negative controls for position 1-29 domain model
# 3x permutate combined_scores, while keeping matrix symmetry
set.seed(1603)
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_1to29.txt"))[Pos1<Pos2,.(Pos1,Pos2,control1 = sample(combined_score),
                                                                          control2 = sample(combined_score),
                                                                          control3 = sample(combined_score))]
PWI2 = switch_double_DT(PWI,list(c("Pos1","Pos2")),c("control1","control2","control3"))
PWI2 = rbind(PWI2,data.table(Pos1=1:29,Pos2=1:29,control1=NA,control2=NA,control3=NA))
write.table(x = PWI2, file = paste0(dataset_dir,"processed_data/DMS_PWI_DC_control.txt"),
            quote = F,row.names = F, col.names = T)
PWI_transformed = deepcontact_transform_basic2d(PWI2,
                                                dataset_dir = dataset_dir,
                                                output_filename = "DMS_PWI_1to29_DC_control_deepcontact.txt",
                                                deepcontact_dir = "where/is/deepcontact/",
                                                prefix = paste0(prefix,"1to29_"))


##################################################################################
##### read distance and secodnary structure info from PDB files (or PSIPRED) #####
##################################################################################

#deposit PDB and PSIPRED files in the dataset subfolders first

## calculate contact maps & secondary structure from PDB file
pairdistances_from_PDB(input_file = paste0(dataset_dir,"dataset/PDB/1k9q.pdb"),
                       dataset_dir = dataset_dir,
                       aa_seq = scan(paste0(dataset_dir,"dataset/WW_sequence.fasta"),what = "character",sep = "\n")[2],
                       idx_pdb_start = 10)

#define betasheets in PDB file 1k9q
beta_sheet_pairing = data.table(pos1_min = c(8,20),
                                pos1_max = c(12,22),
                                pos2_min = c(18,27),
                                pos2_max = c(22,29),
                                type = c("anti-par","anti-par"))
ss_elements = fread(paste0(dataset_dir,"processed_data/PDB_secondary_structure_1k9q_A.txt"))
save(list = c("beta_sheet_pairing","ss_elements"),
     file = paste0(dataset_dir,"processed_data/PDB_secondary_structure_elements_control_1k9q_A.RData"))

# define hbonding in PDB file (1k9q)
# 20F:29T,22N:27T
# 21L:9E,19Y:11A
pdb_beta_sheet = data.table(pos_hn = c(20,29,22,27), pos_o = c(29,20,27,22), strand1 = 2, strand2 = 3, type="anti_par")
pdb_beta_sheet = rbind(pdb_beta_sheet,data.table(pos_hn = c(21,9,19,11), pos_o = c(9,21,11,19), strand1 = 1, strand2 = 2, type="anti_par"))
write.table(x=pdb_beta_sheet,file=paste0(dataset_dir,"processed_data/PDB_beta_sheet_hbonds_1k9q.txt"),quote = F,row.names = F,col.names = T)


#secondary structure from PSIPRED predictions (from PSIPRED webserver v3.3 http://bioinf.cs.ucl.ac.uk/psipred/ )
SS_from_PSIPRED(input_file = paste0(dataset_dir,"dataset/PSIPRED/ww1.ss2"),dataset_dir = dataset_dir)





###########################################
###### evaluate epistasis data (QC) #######
###########################################
### these are some basic scripts to evaluate the dataset and predictions

doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles.txt"))
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))
prefix = "WW_"
contactmap = fread(paste0(dataset_dir,"processed_data/PDB_contactmap_1k9q_A.txt"))

#check how subsets over which positive or negative epistasis can be classified are distributed in single mutant fitness space
epistasis_analytics_subsets_singlemutantspace(doubles,
                                              dataset_dir = dataset_dir,prefix = prefix)

#marginal distribution of # of variants suitable for epistasis classification across position pairs
epistasis_analytics_NumEvars_marginal(doubles,
                                      dataset_dir = dataset_dir,prefix = prefix)

#variants suitable for epistasis classification versus median fitness of single mutants
epistasis_analytics_NumEvars_fitness(doubles,
                                     dataset_dir = dataset_dir,prefix = prefix)

#check the spatial distribution of the number of variants suitable for epistasis classification
epistasis_analytics_NumEvars_spatial(PWI,
                                     dataset_dir = dataset_dir,prefix = prefix)

# cumulative distribution function of epistatic variants as function of distance 
# (if a known contactmap is available)
epistasis_analytics_subsets_CDF(doubles,
                                dataset_dir = dataset_dir,
                                contactmap = contactmap,
                                prefix = prefix,
                                dist_type = "scHAmin", lindist = 5)



###########################################################
######## predict secondary structure from PWI data ########
###########################################################

#### predict secondary structure
predict_secondary_structure_elements(PWI[,.(Pos1,Pos2,
                                            epistasis_score,association_score,combined_score)],
                                     dataset_dir = dataset_dir,
                                     prefix = prefix,
                                     known_SS = paste0(dataset_dir,"processed_data/PDB_secondary_structure_1k9q_A.txt"))

#### predict beta sheets (using PSIPRED predited secondary structure elements as input)
predict_beta_sheets(PWI[,.(Pos1,Pos2,
                           epistasis_score,association_score,combined_score)],
                    input_ss0 = fread(paste0(dataset_dir,"processed_data/PSIPRED_secondary_structure.txt")),
                    dataset_dir = dataset_dir,
                    prefix = paste0(prefix,"PSIPRED_"),
                    known_ss_file = paste0(dataset_dir,"processed_data/PDB_secondary_structure_1k9q_A.txt"),
                    known_bsi_file = paste0(dataset_dir,"processed_data/PDB_beta_sheet_hbonds_1k9q.txt"),
                    debug_this = F,restricted_pairing = T)



##########################################################################################
##### evaluate predicted contacts (top scoring pairs) against reference structure ########
##########################################################################################

#### true positive rate of top contacts + contactmaps + eCDFs
evaluate_contacts_vs_PDB(contacts = PWI[,.(Pos1,Pos2,
                                           posE_enr_mean,negE_enr_mean,epistasis_score,
                                           posE_negE_cor,association_score,combined_score)],
                         contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                         secondary_structure=NA,
                         dataset_dir = dataset_dir,
                         lindist=5,
                         prefix = prefix)

#### minimal number of edges connecting top contacts versus all position pairs
evaluate_contatcs_minimaledges(PWI[,.(Pos1,Pos2,
                                      posE_enr_mean,negE_enr_mean,epistasis_score,
                                      posE_negE_cor,association_score,combined_score)],
                               contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                               dataset_dir = dataset_dir,
                               prefix = prefix,
                               lindist=5,
                               N_contacts = 1,
                               dist_cutoff = 8)

#### interaction scores versus distance in reference structure
score_vs_distance_scatter(contacts = PWI[,.(Pos1,Pos2,
                                            posE_enr_mean,negE_enr_mean,epistasis_score,
                                            posE_negE_cor,association_score,combined_score)],
                          contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                          dataset_dir = dataset_dir,
                          prefix = prefix)




###########################################
####### XPLOR structural modelling ########
###########################################

#use predicted tertiary contacts and secondary structure elements to generate structural models

#### use  DMS data to derive restraints for tertiary contacts & use secondary structure restraints derived from PSIPRED predictions (no beta sheet pairing restraints)
XPLOR_wrapper(input_PWI = merge(PWI[,.(Pos1,Pos2,WT_AA1,WT_AA2,
                                       epistasis_score,association_score,combined_score)],
                                contactmap[,.(Pos1,Pos2,control = scHAmin < 8)],by=c("Pos1","Pos2"),all.x=T),
              SS_mode = "SSonly",
              input_SS_file = "PSIPRED_secondary_structure.txt",
              prefix = "WW_PSIPRED_SSonly_",
              dataset_dir = dataset_dir,
              cores = 15,
              queue = "short-sl7,long-sl7",
              protein_sequence = scan(paste0(dataset_dir,"dataset/WW_sequence.fasta"),sep="\n",what = "character")[2],
              pdb_file = paste0(dataset_dir,"/dataset/PDB/1k9q_model1_mod_manual.pdb"), #this is a modified PDB file to work with XPLOR
              L = c(0.25,0.5,0.75,1,1.25),
              home_dir = "/where/to/temporarily/build/folderstructure/",
              cluster_dir = "cluster/directory/for/XPLOR/",
              login_serveraddress = "mylogin@serveraddress.com",
              debug_this = F) 


########### only for residues 6:29
helper = fread(paste0(dataset_dir,"processed_data/PSIPRED_secondary_structure.txt"))
helper1 = helper[between(Pos,6,29),.(Pos=Pos-5,SS)]
write.table(helper1,file = paste0(dataset_dir,"processed_data/PSIPRED_secondary_structure_6to29.txt"),
            row.names = F,quote = F)
XPLOR_wrapper(input_PWI = merge(PWI[between(Pos1,6,29) & between(Pos2,6,29),
                                    .(Pos1=Pos1-5,Pos2=Pos2-5,
                                      WT_AA1,WT_AA2,
                                      epistasis_score,association_score,combined_score)],
                                contactmap[between(Pos1,6,29) & between(Pos2,6,29),.(Pos1=Pos1-5,Pos2=Pos2-5,control = scHAmin < 8)],
                                by=c("Pos1","Pos2"),all.x=T),
              SS_mode = "SSonly",
              input_SS_file = "PSIPRED_secondary_structure_6to29.txt",
              prefix = "WW_PSIPRED_SSonly_6to29_githubtest",predictorXlength = data.table("combined_score",L=1),
              dataset_dir = dataset_dir,
              cores = 15,
              queue = "short-sl7,long-sl7",
              protein_sequence = paste0(strsplit(scan(paste0(dataset_dir,"dataset/WW_sequence.fasta"),sep="\n",what = "character")[2],"")[[1]][6:29],collapse=""),
              pdb_file = paste0(dataset_dir,"/dataset/PDB/1k9q_model1_6to29.pdb"), #this is a modified PDB file to work with XPLOR
              L = c(0.25,0.5,0.75,1,1.25),
              home_dir = "/where/to/temporarily/build/folderstructure/",
              cluster_dir = "cluster/directory/for/XPLOR/",
              login_serveraddress = "mylogin@serveraddress.com",
              debug_this = F) 




#####################################
##### evaluate XPLOR results ########
#####################################

#this will analyse the "..._variables_results.RData" files generated by the XPLOR simulations
# (copy these back to the dataset directory, e.g. to results/XPLOR/my_prefix_used/, which is the "XPLOR_dir" variable)
# the script will give basic outputs into a "results" sub-directory, such as how models performed across stages, the RMSD/TMscore of the best models etc

analyse_XPLOR_results(XPLOR_dir = paste0(dataset_dir,"results/XPLOR/WW_PSIPRED_SSonly_/"),
                      contactmap = contactmap,
                      eval_L = 0.5,
                      draw_contactmaps = F)

analyse_XPLOR_results(XPLOR_dir = paste0(dataset_dir,"results/XPLOR/WW_PSIPRED_SSonly_6to29_/"),
                      contactmap = contactmap[between(Pos1,6,29) & between(Pos2,6,29),.(Pos1=Pos1-5,Pos2=Pos2-5,scHAmin)],
                      eval_L = 0.5,
                      draw_contactmaps = F)

