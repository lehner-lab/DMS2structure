##############################
##### GB1 from Olson2014 #####
##############################

#this is the pipeline used to analyse the GB1 data from Olson et al. 2014

#first, set the working directory to the DMS2structure folder
setwd("/where/is/DMS2structure/")


#source scripts
filelist = list.files('scripts/')
sapply(paste0('scripts/',filelist),source,.GlobalEnv)

#create the necessary subfolder structure for all results and processed data
dataset_dir = "GB1/"
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

###### load dataset
# this is Olson et al. 2014 Supp. Table S2 from 
# https://www.cell.com/cms/10.1016/j.cub.2014.09.072/attachment/3a36211d-bddd-43e3-bf42-a6721f93a18b/mmc2.xlsx
# table split into three tables for doubles, singles and wildtype; without the "double mutants" header and so on
wildtype = fread(paste0(dataset_dir,"dataset/Olson2014_TableS2_wildtype.txt"), sep = "\t", header = TRUE)
singles = fread(paste0(dataset_dir,"dataset/Olson2014_TableS2_singles.txt"), sep = "\t",header = TRUE)
doubles = fread(paste0(dataset_dir,"dataset/Olson2014_TableS2_doubles.txt"), sep = "\t",header = TRUE)
doubles[,c("V11","V12","V13","V14","V15","V16","V17","V18") := NULL]

## rename coloumns
colnames(wildtype) = c("count_r1_t0","count_r1_t1")
colnames(singles) = c("WT_AA","Pos","Mut","count_r1_t0","count_r1_t1")
colnames(doubles) = c("WT_AA1","Pos1","Mut1","WT_AA2","Pos2","Mut2","count_r1_t0","count_r1_t1","fitness1","fitness2")

###### rearrange doubles such that always Pos1 < Pos2
doubles[Pos1 > Pos2,':=' (Pos1=Pos2,WT_AA1 = WT_AA2,Mut1 = Mut2,fitness1 = fitness2,
                              Pos2=Pos1,WT_AA2 = WT_AA1,Mut2 = Mut1,fitness2 = fitness1)]


###### calculate fitness 
wildtype[,fitness:=0]
singles[,fitness := log(count_r1_t1/count_r1_t0 * (wildtype$count_r1_t0 / wildtype$count_r1_t1))]
#plot wild-type fitness distribution
xd=density(singles$fitness,bw=.15)
plot(xd)
#wild-type peak is not centered at 0; thus try to estimate a correction factor
#fitness peak ~ wildtype peak
correction_factor_wildtype =xd$x[xd$y==max(xd$y)]
#use this as correction factor

###### correct fitness value for this factor
singles[,fitness := NULL]
singles[,fitness := log(count_r1_t1/count_r1_t0 * (wildtype$count_r1_t0 / wildtype$count_r1_t1)) - correction_factor_wildtype]
doubles[,fitness := log(count_r1_t1/count_r1_t0 / (wildtype$count_r1_t1 / wildtype$count_r1_t0)) - correction_factor_wildtype]

# plot fitness distribution
ggplot(rbind(wildtype[,.(fitness,type="wt")],
             singles[,.(fitness,type="singles")],
             doubles[,.(fitness,type="doubles")]),
       aes(fitness,color=type)) +
  geom_density() +
  theme_minimal() +
  labs(x = "log(fitness)")
ggsave(paste0(dataset_dir,"results/preprocessing/fitness_distribution.pdf"),width=4,height = 3)



###### calculate poissonian error
wildtype[,sigma := sqrt(1/count_r1_t1 + 1/count_r1_t0)]
singles[,sigma := sqrt(1/count_r1_t1 + 1/count_r1_t0 + 1/wildtype$count_r1_t1 + 1/wildtype$count_r1_t0)]
doubles[,sigma := sqrt(1/count_r1_t1 + 1/count_r1_t0 + 1/wildtype$count_r1_t1 + 1/wildtype$count_r1_t0)]

# plot fitnes versus error
ggplot(rbind(singles[,.(fitness,sigma,type="singles")],
             doubles[,.(fitness,sigma,type="doubles")]),
       aes(fitness,sigma)) +
  geom_hex() +
  scale_fill_distiller(direction=1,trans="log10") +
  scale_y_log10() +
  facet_wrap( ~ type) +
  theme_minimal() +
  labs(x = "log(fitness)",y = "error of fitness estimate")
ggsave(paste0(dataset_dir,"results/preprocessing/fitness_error_distribution.pdf"),width=8,height=4)



###### pull single mutant fitness and error into double mutant table
doubles[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],.(Pos1,Mut1)]
doubles[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],.(Pos2,Mut2)]
doubles[,sigma1 := singles$sigma[singles$Pos %in% Pos1 & singles$Mut %in% Mut1],by=.(Pos1,Mut1)]
doubles[,sigma2 := singles$sigma[singles$Pos %in% Pos2 & singles$Mut %in% Mut2],by=.(Pos2,Mut2)]



###### mark variants with nonsensical fitness values
wildtype[,is.fitness := fitness > -Inf & !is.na(fitness)]
singles[,is.fitness := fitness > -Inf & !is.na(fitness)]
doubles[,is.fitness := fitness > -Inf & !is.na(fitness)]



###### estimate lower measurement limit of fitness assay 
#from singles via kernel density
xd=density((singles$fitness),bw=.15)
xd$x[xd$y==max(xd$y[xd$x < -4])]
#from doubles that have expected fitness far belowthe lower measurement limit
doubles[is.fitness==T & fitness1 + fitness2 < -8,median(fitness,na.rm=T)]
# nearly identical estimates
lower_bound_F = mean(c(doubles[is.fitness==T & fitness1 + fitness2 < -8,median(fitness,na.rm=T)],
                       xd$x[xd$y==max(xd$y[xd$x < -4])]))



###### define which variants have enough read coverage
# plot read distribution
ggplot(rbind(wildtype[,.(count_r1_t0,count_r1_t1,type="wt")],
             singles[,.(count_r1_t0,count_r1_t1,type="singles")],
             doubles[sample(1:.N,10000),.(count_r1_t0,count_r1_t1,type="doubles")]),
       aes(count_r1_t0,count_r1_t1,color=type)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(x = "read counts in input library",
       y = "read counts in output library")
ggsave(paste0(dataset_dir,"results/preprocessing/readcounts_input_output.pdf"),width=6,height=6)
# > double mutants with low input coverage run into limitations

wildtype[,is.reads0 := TRUE]
singles[,is.reads0 := TRUE]
# only throw away variants with zero output counts if unclear where above lower_fitness_bound their fitness would be
# only applies to doubles
lower_read_cut = 200
doubles[,is.reads0 := count_r1_t0 > 10 & (count_r1_t1 >= 1 | count_r1_t0 >= lower_read_cut)]
# plot this
ggplot() +
  geom_hex(data=doubles[is.fitness==T],aes(x=count_r1_t0,y=fitness),bins=75) +
  scale_fill_gradient(low="gray95", high = "dodgerblue4") +
  scale_x_log10(breaks=10^seq(0,5,1)) +
  geom_hline(yintercept = lower_bound_F,color="black",linetype=2) +
  #plot input read restrictions
  geom_line(data=data.frame(x=c(10.5,10.5),y=c(log(0.75/10.5*wildtype[,count_r1_t0/count_r1_t1]),2.5)), aes(x,y),color="red") +
  geom_line(data=data.frame(x=c(10.5,lower_read_cut),
                            y=c(log(0.75/10.5*wildtype[,count_r1_t0/count_r1_t1]),
                                log(0.75/lower_read_cut*wildtype[,count_r1_t0/count_r1_t1]))), aes(x,y),color="red") +
  geom_line(data=data.frame(x=c(lower_read_cut,lower_read_cut),
                            y=c(log(0.75/lower_read_cut*wildtype[,count_r1_t0/count_r1_t1]),-8)), aes(x,y),color="red") +
  geom_text(data = data.frame(x=11,y=2,
                              label = paste0(">> ",round(doubles[is.reads0==T,.N]/doubles[,.N],digits=2)*100,"%")),
            aes(x,y,label=label),hjust=0) +
  coord_cartesian(xlim=c(1,10^5),expand=c(0,0)) +
  scale_y_continuous(breaks=seq(-10,2.5,2.5)) +
  theme_classic() +
  labs(y="fitness (log)", x="sequencing reads input",fill = "# variants")
ggsave(paste0(dataset_dir,"results/preprocessing/readcounts_input_versus_fitness.pdf"),width=6,height=4)



###### reorder doubles table
doubles = doubles[,.SD,,.SDcols = c("Pos1","Pos2","Mut1","Mut2","WT_AA1","WT_AA2",
                      "count_r1_t0","count_r1_t1","is.fitness","is.reads0",
                      "fitness1","fitness2","sigma1","sigma2",
                      "fitness","sigma")]




###### write data tables to processed_data folder
write.table(x = wildtype, file = paste0(dataset_dir,"processed_data/DMS_wildtype.txt"),
            quote = F,row.names = F, col.names = T)
write.table(x = singles, file = paste0(dataset_dir,"processed_data/DMS_singles.txt"),
            quote = F,row.names = F, col.names = T)
write.table(x = doubles, file = paste0(dataset_dir,"processed_data/DMS_doubles_preE.txt"),
            quote = F,row.names = F, col.names = T)



##################################################################
### calculate epistasis null model and call pos./neg.epistasis ###
##################################################################

#call epistatic interactions
# doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles_preE.txt"))
doubles = call_epistasis_binary(doubles,
                                lower_bound_F,
                                dataset_dir = dataset_dir,
                                prefix = "GB1_")


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

# PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))
PWI_transformed = deepcontact_transform_basic2d(PWI[,.(Pos1,Pos2,
                                              epistasis_score,association_score,combined_score,
                                              combined_score_negE,combined_score_posE)],
                                        dataset_dir = dataset_dir,
                                        deepcontact_dir = "where/is/deepcontact/",
                                        prefix = "GB1_")

### negative controls for DeepContact
# 3x permutate combined_scores, while keeping matrix symmetry
set.seed(1603)
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))[Pos1<Pos2,.(Pos1,Pos2,control1 = sample(combined_score),
                                                                          control2 = sample(combined_score),
                                                                          control3 = sample(combined_score))]
PWI2 = switch_double_DT(PWI,list(c("Pos1","Pos2")),c("control1","control2","control3"))
PWI2 = rbind(PWI2,data.table(Pos1=2:56,Pos2=2:56,control1=NA,control2=NA,control3=NA))
write.table(x = PWI2, file = paste0(dataset_dir,"processed_data/DMS_PWI_DC_control.txt"),
            quote = F,row.names = F, col.names = T)
PWI_transformed = deepcontact_transform_basic2d(PWI2,
                                           dataset_dir = dataset_dir,
                                           output_filename = "DMS_PWI_DC_control_deepcontact.txt",
                                           deepcontact_dir = "where/is/deepcontact/",
                                           prefix = "GB1_")





##################################################################################
##### read distance and secodnary structure info from PDB files (or PSIPRED) #####
##################################################################################

#deposit PDB and PSIPRED files in the dataset subfolders first

## calculate contact maps & secondary structure from PDB file
pairdistances_from_PDB(input_file = paste0(dataset_dir,"dataset/PDB/1pga.pdb"),
                       dataset_dir = dataset_dir,
                       aa_seq = scan(paste0(dataset_dir,"dataset/GB1_sequence.fasta"),what = "character")[2])

## define betasheets in PDB file 1pga (for use as control in structural simulations) & combine with "true" secondary structure elements from PDB file
#this is extracted from PyMOL
beta_sheet_pairing = data.table(pos1_min = c(1,3,42),
                                pos1_max = c(9,9,46),
                                pos2_min = c(12,50,51),
                                pos2_max = c(20,56,55),
                                type = c("anti-par","par","anti-par"))
ss_elements = fread(paste0(dataset_dir,"processed_data/PDB_secondary_structure_1pga_A.txt"))
save(list = c("beta_sheet_pairing","ss_elements"),
     file = paste0(dataset_dir,"processed_data/PDB_secondary_structure_elements_control_1pga_A.RData"))

## define hbonding in PDB file (1pga)
pdb_beta_sheet = data.table(pos_hn = c(20,3,18,5,16,7,14,9), pos_o = c(1,18,3,16,5,14,7,12), strand1 = 1, strand2 = 2, type="anti_par")
pdb_beta_sheet = rbind(pdb_beta_sheet,data.table(pos_hn = c(4,52,6,54,8,56), pos_o = c(50,4,52,6,54,8), strand1 = 1, strand2=4, type="par"))
pdb_beta_sheet = rbind(pdb_beta_sheet,data.table(pos_hn = c(42,55,44,53,46,51), pos_o = c(55,42,53,44,51,46), strand1 = 3, strand2=4, type="anti_par"))
write.table(x=pdb_beta_sheet,file=paste0(dataset_dir,"processed_data/PDB_beta_sheet_hbonds_1pga.txt"),quote = F,row.names = F,col.names = T)

#secondary structure from PSIPRED predictions (from PSIPRED webserver v3.3 http://bioinf.cs.ucl.ac.uk/psipred/ )
SS_from_PSIPRED(input_file = paste0(dataset_dir,"dataset/PSIPRED/gb1.ss2"),dataset_dir = dataset_dir)



###########################################
###### evaluate epistasis data (QC) #######
###########################################
### these are some basic scripts to evaluate the dataset and predictions

doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles.txt"))
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))
prefix = "GB1_"
contactmap = fread(paste0(dataset_dir,"processed_data/PDB_contactmap_1pga_A.txt"))

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
                                  known_SS = paste0(dataset_dir,"processed_data/PDB_secondary_structure_1pga_A.txt"))
  
#### predict beta sheets
predict_beta_sheets(PWI[,.(Pos1,Pos2,
                           epistasis_score,association_score,combined_score)],
                    input_ss0 = fread(paste0(dataset_dir,"processed_data/",prefix,"secondary_structure_prediction.txt")),
                    dataset_dir = dataset_dir,
                    prefix = prefix,
                    known_ss_file = paste0(dataset_dir,"processed_data/PDB_secondary_structure_1pga_A.txt"),
                    known_bsi_file = paste0(dataset_dir,"processed_data/PDB_beta_sheet_hbonds_1pga.txt"))





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

#these scripts are multi-level wrapper functions to run XPLOR-NIH on a compute cluster
#they locally create/modify the necessary scripts and folder structure, and scp this to the cluster, then remotely execute parallel jobs
#this might take quite some work to adapt to other working environments


#### use only DMS data to derive restraints for tertiary contacts & secondary structure elements (including beta sheet pairing)
XPLOR_wrapper(input_PWI =  merge(PWI[,.(Pos1,Pos2,WT_AA1,WT_AA2,
                                        epistasis_score,association_score,combined_score)],
                                 contactmap[,.(Pos1,Pos2,control = scHAmin < 8)], #as positive/negative controls
                                 by=c("Pos1","Pos2"),all.x=T),
              SS_mode = "SSsheets",
              input_SS_file = "GB1_",
              prefix = paste0(prefix,"SSsheets_"),
              dataset_dir = dataset_dir,
              cores = 15,
              queue = "short-sl7,long-sl7",
              protein_sequence = scan(paste0(dataset_dir,"dataset/GB1_sequence.fasta"),what = "character")[2],
              pdb_file = paste0(dataset_dir,"dataset/PDB/g_xray.pdb"),
              home_dir = "/where/to/temporarily/build/folderstructure/",
              cluster_dir = "cluster/directory/for/XPLOR/",
              login_serveraddress = "mylogin@serveraddress.com",
              debug_this = F)

#### use  DMS data to derive restraints for tertiary contacts & use secondary structure restraints derived from PSIPRED predictions (no beta sheet pairing restraints)
# +compare to deepcontact transformed scores
helper1 = merge(PWI[,.(Pos1,Pos2,WT_AA1,WT_AA2,
             epistasis_score,association_score,combined_score,
             combined_score_negE,combined_score_posE)],
      contactmap[,.(Pos1,Pos2,control = scHAmin < 8)],by=c("Pos1","Pos2"),all.x=T)
PWI_DC = fread(paste0(dataset_dir,"processed_data/DMS_PWI_deepcontact.txt"))
helper2 = merge(helper1,PWI_DC[,.(Pos1,Pos2,
                                  epistasis_score_DC = epistasis_score,
                                  association_score_DC = association_score,
                                  combined_score_DC = combined_score,
                                  combined_score_negE_DC = combined_score_negE,
                                  combined_score_posE_DC = combined_score_posE)],
                by=c("Pos1","Pos2"),all=T)

XPLOR_wrapper(input_PWI = helper2,
              SS_mode = "SSonly",
              input_SS_file = "PSIPRED_secondary_structure.txt",
              prefix = paste0(prefix,"PSIPRED_SSonly_"),
              dataset_dir = dataset_dir,
              cores = 15,
              queue = "short-sl7,long-sl7",
              protein_sequence = scan(paste0(dataset_dir,"dataset/GB1_sequence.fasta"),what = "character")[2],
              pdb_file = paste0(dataset_dir,"dataset/PDB/g_xray.pdb"),
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

analyse_XPLOR_results(XPLOR_dir = paste0(dataset_dir,"results/XPLOR/GB1_SSsheets_/"),
                      contactmap = fread(paste0(dataset_dir,"processed_data/PDB_contactmap_1pga_A.txt")),
                      draw_contactmaps = F)


analyse_XPLOR_results(XPLOR_dir = paste0(dataset_dir,"results/XPLOR/GB1_PSIPRED_SSonly_/"),
                      contactmap = fread(paste0(dataset_dir,"processed_data/PDB_contactmap_1pga_A.txt")),
                      draw_contactmaps = F)

