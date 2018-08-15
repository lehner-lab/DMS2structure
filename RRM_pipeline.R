#######################################
##### RRM domain from Melamed2013 ##### of PAB1 yeast gene
#######################################

#this is the pipeline used to analyse the RRM domain data from Melamed et al. 2013

#first, set the working directory to the DMS2structure folder
setwd("/where/is/DMS2structure/")


#source scripts
filelist = list.files('scripts/')
sapply(paste0('scripts/',filelist),source,.GlobalEnv)

#create the necessary subfolder structure for all results and processed data
dataset_dir = "RRM/"
create_directory_structure(dataset_dir)
#then save this script in the dataset_dir
#and paste all necessary source data into dataset_dir/dataset

#load required packages
require(data.table)
require(ggplot2)
require(cowplot)

#############################################################################################
##### preprocess data (calculate fitness scores and errors, set quality thresholds etc) #####
#############################################################################################

################# read data from Supplementary Table 5 from Melamed2013
dataset = fread(paste0(dataset_dir,"dataset/Supplementary_Table_5_doubles.txt"), sep = "\t",header = TRUE)

# extract position and amino acids
dataset[,Pos1:=sapply(strsplit(dataset[,as.character(seqID_X)],"-"),FUN = function(X){as.integer(X[1])})]
dataset[,Mut1:=sapply(strsplit(dataset[,as.character(seqID_X)],"-"),FUN = function(X){X[2]})]

dataset[,Pos2:=sapply(strsplit(dataset[,as.character(seqID_Y)],"-"),FUN = function(X){as.integer(X[1])})]
dataset[,Mut2:=sapply(strsplit(dataset[,as.character(seqID_Y)],"-"),FUN = function(X){X[2]})]

# indicate library
dataset[Pos1<=150,section := 1]
dataset[between(Pos1,151,175),section := 2]
dataset[Pos1>175,section := 3]


#investigate dependency on input counts
ggplot(dataset,aes(Input_reads,log(XY_Enrichment_score))) + 
  geom_hex() +
  scale_x_log10() +
  scale_fill_gradient(trans="log",breaks=c(1,10,100)) +
  facet_wrap(~section)
ggsave(paste0(dataset_dir,"/results/preprocessing/RRM_Inputreads_fitness.pdf"),width=8,height=3)

### should treat libraries separately, develop independent error estimates
### output count should be Input_reads * XY_enrichment_score
#define abbrev. variables
dataset = dataset[,.(Pos1,Pos2,Mut1,Mut2,section,dist = Physical_distance,
                         count_e1_s0 = Input_reads,count_e1_s1 = Input_reads * X_Enrichment_score,
                         fitness1=log(X_Enrichment_score),fitness2=log(Y_Enrichment_score),fitness=log(XY_Enrichment_score))]

# calculate Poissonian error
dataset[,sigma := sqrt(1/count_e1_s0 + 1/count_e1_s1)]

ggplot(dataset,aes(fitness,sigma)) + 
  geom_hex() + 
  scale_y_log10() + 
  scale_fill_continuous(trans = "log10") +
  facet_grid( ~ section)
ggsave(paste0(dataset_dir,"/results/preprocessing/RRM_fitness_sigma.pdf"),width=8,height=3)


# define STOP variants
dataset[,STOP := Mut1 == "*" | Mut2 == "*"]

# plot fitness distribution + STOPs
ggplot(dataset,aes(fitness,color=factor(section),linetype=STOP)) + 
  geom_density()
ggsave(paste0(dataset_dir,"/results/preprocessing/RRM_fitness_distribution.pdf"),width=5,height=3)
# >> all libraries are different

#estimate lower fitness bounds from STOP variants
dataset[,lower_bound_F := weighted.mean(x = dataset[section == unique(unlist(.SD)) & STOP == TRUE,fitness],
                                        w = dataset[section == unique(unlist(.SD)) & STOP == TRUE,1/sigma^2]),
        section,.SDcols = "section"]
unique(dataset[,.(lower_bound_F,section)][order(section)])

## mark STOP variants and variants with enough reads (all)
dataset[,is.fitness := !STOP]
dataset[,is.reads0 := T]

#give dummy sigmas for single mutants (no read counts for singles available)
dataset[,sigma1 := 0.01]
dataset[,sigma2 := 0.01]

#correct Positions
dataset[,Pos1 := Pos1 - 125]
dataset[,Pos2 := Pos2 - 125]

#add WT_AA
WT_aaseq = scan(paste0(dataset_dir,"dataset/RRM_domain_sequence.fasta"),what="character",sep="\n")[2]
WT_aaseq_split = data.table(Pos = 1:75,WT_AA = strsplit(WT_aaseq,"")[[1]])

dataset[,WT_AA1 := WT_aaseq_split[Pos == Pos1,WT_AA],Pos1]
dataset[,WT_AA2 := WT_aaseq_split[Pos == Pos2,WT_AA],Pos2]


# reorder double table
doubles = dataset[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,
                     section,count_e1_s0,count_e1_s1,STOP,is.fitness,is.reads0,
                     fitness,sigma,fitness1,sigma1,fitness2,sigma2,lower_bound_F)]

#### save data.tables
write.table(x = doubles, file = paste0(dataset_dir,"processed_data/DMS_doubles_preE.txt"),
            quote = F,row.names = F, col.names = T)


### also save extracted pair-distances from dataset table
ggplot(unique(dataset[,.(Pos1,Pos2,dist,section)]),aes(Pos1,Pos2,fill=dist<8)) +
  geom_raster() +
  scale_fill_manual(values=c("grey95","grey25"),na.value = "red")
ggsave(paste0(dataset_dir,"/results/preprocessing/RRM_contactmap_Melamed2013.pdf"),width=6,height=5)

#only section 2 has a decent number of off-diagonal contacts
# !!! this is missing position 14 of the second section (position 39 in absolute terms)

write.table(x = unique(dataset[section == 1,.(Pos1,Pos2,WT_AA1,WT_AA2,dist,section)]), file = paste0(dataset_dir,"processed_data/contactmap_RRM_sec1.txt"),
            quote = F,row.names = F, col.names = T)

write.table(x = unique(dataset[section == 2 & !is.na(dist),.(Pos1 = Pos1-25,Pos2 = Pos2-25,WT_AA1,WT_AA2,dist,section)]), file = paste0(dataset_dir,"processed_data/contactmap_RRM_sec2.txt"),
            quote = F,row.names = F, col.names = T)

write.table(x = unique(dataset[section == 3,.(Pos1 = Pos1-50,Pos2 = Pos2-50,WT_AA1,WT_AA2,dist,section)]), file = paste0(dataset_dir,"processed_data/contactmap_RRM_sec3.txt"),
            quote = F,row.names = F, col.names = T)



##################################################################
### calculate epistasis null model and call pos./neg.epistasis ###
##################################################################

## given that only section 2 has off-diagonal contacts, focus on this
doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles_preE.txt"))
doubles2 = copy(doubles[section==2])
doubles2[,':=' (Pos1 = Pos1-25,Pos2 = Pos2-25)]
doubles = call_epistasis_binary(doubles2,
                                unique(doubles2$lower_bound_F),
                                dataset_dir = dataset_dir,
                                prefix = "RRM_sec2_",
                                output_filename = "DMS_doubles_sec2.txt",
                                epistasis_error_from_slopes = F)

#############################################
### calculate pairwise interaction scores ###
#############################################

doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles_sec2.txt"))
PWI = calculate_pairwise_interaction_scores(doubles,
                                            dataset_dir = dataset_dir,
                                            output_filename = "DMS_PWI_sec2.txt",
                                            detailed = F)


#########################################
### deep contact transform PWI scores ###
#########################################

PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_sec2.txt"))
PWI_transformed = deepcontact_transform_basic2d(PWI,
                                                dataset_dir = dataset_dir,
                                                deepcontact_dir = "where/is/deepcontact/",
                                                output_filename = "DMS_PWI_sec2_deepcontact.txt",
                                                prefix = "RRM_sec2_")

### negative controls for DeepContact
# 3x permutate combined_scores, while keeping matrix symmetry
set.seed(1603)
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_sec2.txt"))[Pos1<Pos2,.(Pos1,Pos2,control1 = sample(combined_score),
                                                                          control2 = sample(combined_score),
                                                                          control3 = sample(combined_score))]
PWI2 = switch_double_DT(PWI,list(c("Pos1","Pos2")),c("control1","control2","control3"))
PWI2 = rbind(PWI2,data.table(Pos1=1:25,Pos2=1:25,control1=NA,control2=NA,control3=NA))
write.table(x = PWI2, file = paste0(dataset_dir,"processed_data/DMS_PWI_sec2_DC_control.txt"),
            quote = F,row.names = F, col.names = T)
PWI_transformed = deepcontact_transform_basic2d(PWI2,
                                                dataset_dir = dataset_dir,
                                                output_filename = "DMS_PWI_sec2_DC_control_deepcontact.txt",
                                                deepcontact_dir = "where/is/deepcontact/",
                                                prefix = "RRM_sec2_")

##################################################################################
##### read distance and secodnary structure info from PDB files (or PSIPRED) #####
##################################################################################

#deposit PDB and PSIPRED files in the dataset subfolders first

## calculate contact maps & secondary structure from PDB file
pairdistances_from_PDB(input_file = paste0(dataset_dir,"dataset/PDB/1cvj.pdb"),
                       dataset_dir = dataset_dir,
                       aa_seq = paste0(WT_aaseq_split[Pos %in% c(26:38,40:50),WT_AA],collapse=""), ### as in Melamed2013, omit position 39
                       idx_pdb_start = 124,idx_DMS_start = 1,
                       given_chainids = "A",suffix = "_sec2",debug_this = F)

#define betasheets in PDB file 1cvj for secction 2
beta_sheet_pairing = data.table(pos1_min = c(126-123),
                                pos1_max = c(130-123),
                                pos2_min = c(141-123),
                                pos2_max = c(145-123),
                                type = c("anti-par"))
ss_elements = fread(paste0(dataset_dir,"processed_data/PDB_secondary_structure_1cvj_A_sec2.txt"))
save(list = c("beta_sheet_pairing","ss_elements"),
     file = paste0(dataset_dir,"processed_data/PDB_secondary_structure_elements_control_1cvj_A_sec2.RData"))

#define hbonding in PDB file (1cvj, section 2)
pdb_beta_sheet = data.table(pos_hn = c(127-123,144-123,129-123,142-123), 
                            pos_o = c(144-123,127-123,142-123,129-123), strand1 = 1, strand2 = 2, type="anti_par")
write.table(x=pdb_beta_sheet,file=paste0(dataset_dir,"processed_data/PDB_beta_sheet_hbonds_1cvj_A_sec2.txt"),quote = F,row.names = F,col.names = T)

SS_from_PSIPRED(input_file = paste0(dataset_dir,"dataset/PSIPRED/PAB1.ss2"),dataset_dir = dataset_dir)


######################################################
### correct sec2 data to fit with 1cvj contactmap  ###
######################################################

PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_sec2.txt"))
PWI_corr = PWI[Pos1 != 14 & Pos2 != 14] #remove position 14 for which there is no distance information
PWI_corr[,':=' (Pos1 = ifelse(Pos1 < 14,Pos1,Pos1-1),Pos2 = ifelse(Pos2 < 14,Pos2,Pos2-1))]
write.table(x = PWI_corr, file = paste0(dataset_dir,"processed_data/DMS_PWI_sec2_corr14.txt"),
            quote = F,row.names = F, col.names = T)

PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_sec2_deepcontact.txt"))
PWI_corr = PWI[Pos1 != 14 & Pos2 != 14]
PWI_corr[,':=' (Pos1 = ifelse(Pos1 < 14,Pos1,Pos1-1),Pos2 = ifelse(Pos2 < 14,Pos2,Pos2-1))]
write.table(x = PWI_corr, file = paste0(dataset_dir,"processed_data/DMS_PWI_sec2_corr14_deepcontact.txt"),
            quote = F,row.names = F, col.names = T)


### for randomized DC controls
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_sec2_DC_control.txt"))
PWI_corr = PWI[Pos1 != 14 & Pos2 != 14]
PWI_corr[,':=' (Pos1 = ifelse(Pos1 < 14,Pos1,Pos1-1),Pos2 = ifelse(Pos2 < 14,Pos2,Pos2-1))]
write.table(x = PWI_corr, file = paste0(dataset_dir,"processed_data/DMS_PWI_sec2_corr14_DC_control.txt"),
            quote = F,row.names = F, col.names = T)

PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_sec2_DC_control_deepcontact.txt"))
PWI_corr = PWI[Pos1 != 14 & Pos2 != 14]
PWI_corr[,':=' (Pos1 = ifelse(Pos1 < 14,Pos1,Pos1-1),Pos2 = ifelse(Pos2 < 14,Pos2,Pos2-1))]
write.table(x = PWI_corr, file = paste0(dataset_dir,"processed_data/DMS_PWI_sec2_corr14_DC_control_deepcontact.txt"),
            quote = F,row.names = F, col.names = T)




###########################################
###### evaluate epistasis data (QC) #######
###########################################
### these are some basic scripts to evaluate the dataset and predictions

doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles_sec2.txt"))
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_sec2_corr14.txt"))
prefix = "RRM_sec2_"
contactmap = fread(paste0(dataset_dir,"processed_data/PDB_contactmap_1cvj_A_sec2.txt"))

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
predict_secondary_structure_elements(PWI,
                                     dataset_dir = dataset_dir,
                                     prefix = prefix,
                                     known_SS = paste0(dataset_dir,"processed_data/PDB_secondary_structure_1cvj_A_sec2.txt"))

#### predict beta sheets given PSIPRED secondary structure prediction
#first modify PSIPRED input 
#i.e. isolate section 2 positions (26-50), then remove position 14 of this section
psipred = fread(paste0(dataset_dir,"processed_data/PSIPRED_secondary_structure.txt"))[between(Pos,26,50),.(Pos=Pos-25,SS)][Pos!=14,.(Pos = ifelse(Pos<14,Pos,Pos-1),SS)]
predict_beta_sheets(PWI,
                    input_ss0 = psipred,
                    dataset_dir = dataset_dir,
                    prefix = prefix,
                    known_ss_file = paste0(dataset_dir,"processed_data/PDB_secondary_structure_1cvj_A_sec2.txt"),
                    known_bsi_file = paste0(dataset_dir,"processed_data/PDB_beta_sheet_hbonds_1cvj_A_sec2.txt"),
                    restricted_pairing = T)



##########################################################################################
##### evaluate predicted contacts (top scoring pairs) against reference structure ########
##########################################################################################

#### true positive rate of top contacts + contactmaps + eCDFs
evaluate_contacts_vs_PDB(contacts = PWI,
                         contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                         secondary_structure=NA,
                         dataset_dir = dataset_dir,
                         lindist=5,
                         prefix = prefix)

#### minimal number of edges connecting top contacts versus all position pairs
evaluate_contatcs_minimaledges(PWI,
                               contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                               dataset_dir = dataset_dir,
                               prefix = prefix,
                               lindist=5,
                               N_contacts = 1,
                               dist_cutoff = 8)

#### interaction scores versus distance in reference structure
score_vs_distance_scatter(contacts = PWI,
                          contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                          dataset_dir = dataset_dir,
                          prefix = prefix)


################################
### contacts vs PDB ############
################################

evaluate_contacts_vs_PDB(contacts = PWI[,.(Pos1,Pos2,
                                           epistasis = posE_negE_enr,association = posE_negE_pcor,combined = posE_negE_pcor_enr,
                                           negE_epistasis = negE_enr,negE_association = negE_pcor, negE_combined = negE_pcor_enr,
                                           posE_epistasis = posE_enr,posE_association = posE_pcor, posE_combined = posE_pcor_enr)],
                         contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                         secondary_structure=NA,
                         dataset_dir = dataset_dir,lindist=5,
                         prefix = prefix)
evaluate_contacts_vs_PDB(contacts = PWI[Pos1 != 14 & Pos2 != 14,.(Pos1 = ifelse(Pos1 < 14,Pos1,Pos1-1),
                                                                  Pos2 = ifelse(Pos2 < 14,Pos2,Pos2-1),
                                           epistasis = posE_negE_enr,association = posE_negE_pcor,combined = posE_negE_pcor_enr,
                                           negE_epistasis = negE_enr,negE_association = negE_pcor, negE_combined = negE_pcor_enr,
                                           posE_epistasis = posE_enr,posE_association = posE_pcor, posE_combined = posE_pcor_enr)],
                         contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                         secondary_structure=NA,
                         dataset_dir = dataset_dir,lindist=5,
                         prefix = paste0(prefix,"m14_"))

contacts_vs_distance_scatter(contacts = PWI[,.(Pos1,Pos2,
                                               epistasis = posE_negE_enr,association = posE_negE_pcor,combined = posE_negE_pcor_enr,
                                               negE_epistasis = negE_enr,negE_association = negE_pcor, negE_combined = negE_pcor_enr,
                                               posE_epistasis = posE_enr,posE_association = posE_pcor, posE_combined = posE_pcor_enr)],
                             contactmap = contactmap[,.(Pos1,Pos2,dist)],
                             dataset_dir = dataset_dir,
                             prefix = prefix)

