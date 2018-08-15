############################
##### GB1 downsampling #####
############################

#this is the pipeline used to analyse the downsampled versions of the GB1 data from Olson et al. 2014

#first, set the working directory to the DMS2structure folder
setwd("/where/is/DMS2structure/")

#source scripts
filelist = list.files('scripts/')
sapply(paste0('scripts/',filelist),source,.GlobalEnv)

#create the necessary subfolder structure for all results and processed data
dataset_dir = "GB1_downsampling/"
create_directory_structure(dataset_dir)
#then save this script in the dataset_dir
#and paste all necessary source data into dataset_dir/dataset

#load required packages
require(data.table)
require(ggplot2)
require(cowplot)
require(GGally)
theme_set(theme_minimal())
require(seqinr)


###############################
##### READ downsampling #######
###############################

read_downsampling = c(0.25,0.1,0.025)

for (RD in seq_along(read_downsampling)) {
  
  #load original data
  wildtype = fread("GB1/dataset/Olson2014_TableS2_wildtype.txt", sep = "\t", header = TRUE)
  singles = fread("GB1/dataset/Olson2014_TableS2_singles.txt", sep = "\t",header = TRUE)
  doubles = fread("GB1/dataset/Olson2014_TableS2_doubles.txt", sep = "\t",header = TRUE)
  doubles[,c("V11","V12","V13","V14","V15","V16","V17","V18") := NULL]
  
  # rename coloumns
  colnames(wildtype) = c("count_r1_t0","count_r1_t1")
  colnames(singles) = c("WT_AA","Pos","Mut","count_r1_t0","count_r1_t1")
  colnames(doubles) = c("WT_AA1","Pos1","Mut1","WT_AA2","Pos2","Mut2","count_r1_t0","count_r1_t1","fitness1","fitness2")
  
  
  ## rearrange doubles$GB1 such that always Pos1 < Pos2
  doubles[Pos1 > Pos2,':=' (Pos1=Pos2,WT_AA1 = WT_AA2,Mut1 = Mut2,fitness1 = fitness2,
                            Pos2=Pos1,WT_AA2 = WT_AA1,Mut2 = Mut1,fitness2 = fitness1)]
  
  
  ########################
  ##### downsample #######
  ########################
  set.seed(1603)
  wildtype[,count_r1_t0 := rbinom(1,size = count_r1_t0,prob = read_downsampling[RD])]
  wildtype[,count_r1_t1 := rbinom(1,size = count_r1_t1,prob = read_downsampling[RD])]
  
  singles[,count_r1_t0 := rbinom(1,size = count_r1_t0,prob = read_downsampling[RD]),.(Pos,Mut)]
  singles[,count_r1_t1 := rbinom(1,size = count_r1_t1,prob = read_downsampling[RD]),.(Pos,Mut)]
  
  doubles[,count_r1_t0 := rbinom(1,size = count_r1_t0,prob = read_downsampling[RD]),.(Pos1,Pos2,Mut1,Mut2)]
  doubles[,count_r1_t1 := rbinom(1,size = count_r1_t1,prob = read_downsampling[RD]),.(Pos1,Pos2,Mut1,Mut2)]
  
  
  ## calculate fitness
  wildtype[,fitness:=0]
  singles[,fitness := log(count_r1_t1/count_r1_t0 * (wildtype$count_r1_t0 / wildtype$count_r1_t1))]
  
  #wild-type correction factor
  xd=density(singles$fitness,bw=.15)
  # plot(xd)
  #fitness peak ~ wildtype peak
  ## both fits give similar result for upper mode, exp(0.193) and exp(0.169)
  correction_factor_wildtype =xd$x[xd$y==max(xd$y)]
  
  #correct fitness value for this factor
  singles[,fitness := NULL]
  singles[,fitness := log(count_r1_t1/count_r1_t0 * (wildtype$count_r1_t0 / wildtype$count_r1_t1)) - correction_factor_wildtype]
  doubles[,fitness := log(count_r1_t1/count_r1_t0 / (wildtype$count_r1_t1 / wildtype$count_r1_t0)) - correction_factor_wildtype]
  
  # calculate standard-error of fitness values given read counts
  wildtype[,sigma := sqrt(1/count_r1_t1 + 1/count_r1_t0)]
  singles[,sigma := sqrt(1/count_r1_t1 + 1/count_r1_t0 + 1/wildtype$count_r1_t1 + 1/wildtype$count_r1_t0)]
  doubles[,sigma := sqrt(1/count_r1_t1 + 1/count_r1_t0 + 1/wildtype$count_r1_t1 + 1/wildtype$count_r1_t0)]
  
  # transfer single fitness/error values to doubles data.table
  doubles[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],.(Pos1,Mut1)]
  doubles[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],.(Pos2,Mut2)]
  doubles[,sigma1 := singles$sigma[singles$Pos %in% Pos1 & singles$Mut %in% Mut1],by=.(Pos1,Mut1)]
  doubles[,sigma2 := singles$sigma[singles$Pos %in% Pos2 & singles$Mut %in% Mut2],by=.(Pos2,Mut2)]
  
  #mark variants with nonsensical fitness values
  wildtype[,is.fitness := fitness > -Inf & !is.na(fitness)]
  singles[,is.fitness := fitness > -Inf & !is.na(fitness)]
  doubles[,is.fitness := fitness > -Inf & !is.na(fitness) & fitness1 > -Inf & fitness2 > -Inf]
  
  # define which variants have enough reads
  wildtype[,is.reads0 := TRUE]
  singles[,is.reads0 := TRUE]
  # only throw away variants with zero output counts if unclear where above lower_fitness_bound their fitness would be
  # only applies to doubles
  lower_read_cut = 200
  doubles[,is.reads0 := count_r1_t0 > 10 & (count_r1_t1 >= 1 | count_r1_t0 >= lower_read_cut)]
  
  # rearrange doubles data.table
  doubles = doubles[,.SD,,.SDcols = c("Pos1","Pos2","Mut1","Mut2","WT_AA1","WT_AA2",
                                      "count_r1_t0","count_r1_t1","is.fitness","is.reads0",
                                      "fitness1","fitness2","sigma1","sigma2",
                                      "fitness","sigma")]
  
  # save doubles data.table
  write.table(x = doubles, file = paste0(dataset_dir,"processed_data/DMS_doubles_preE_RD",read_downsampling[RD],".txt"),
              quote = F,row.names = F, col.names = T)
}



############################
##### doped dataset  #######
############################
# only allow AA mutations 1nt hamming distance away from coding sequence
# load a codon<>AA conversion table
AA_codon_conversion = read.table(paste0(dataset_dir,"dataset/amino_acid_codon_conversion.txt"),sep = "\t",header = T)
# load the coding sequence of G protein B1 domain
Gprot_nuc_seq = read.fasta(paste0(dataset_dir,"dataset/GB1_CDS_nt"))
## 229:282 corresponds to 3:56 in GB1_Olson AA sequence
# seqinr::translate(Gprot_nuc_seq$`M13825.1:578-1924`)[227:282]
# seqinr::translate(Gprot_nuc_seq$`M13825.1:578-1924`[(227*3-2):(282*3)])
acgt = c("a","c","g","t")
GB1_nuc_seq = Gprot_nuc_seq$`M13825.1:578-1924`[(227*3-2):(282*3)]

#fix position 2 to Q/Glutamine, Codon CAA
GB1_nuc_seq[4:6] = c("c","a","a")
GB1_doped_mut = data.table(Pos = rep(0,55*12), WTaa = rep("",55*12), Mut = rep("",55*12))
for (aa_pos in 2:56) {
  for (nt_pos in 1:3) {
    for (nt in 1:4) {
      GB1_doped_mut[(aa_pos-2)*12 + (nt_pos-1)*4 + nt,Pos:= aa_pos]
      wt = GB1_nuc_seq[((aa_pos-1)*3+1) : ((aa_pos-1)*3+3)]
      GB1_doped_mut[(aa_pos-2)*12 + (nt_pos-1)*4 + nt,WTaa := seqinr::translate(wt)]
      mutated = wt
      mutated[nt_pos] = acgt[nt]
      GB1_doped_mut[(aa_pos-2)*12 + (nt_pos-1)*4 + nt,Mut := seqinr::translate(mutated)]
    } 
  }
}
print(paste("possible NT1 mutations:",nrow(GB1_doped_mut)))
#get rid of duplicate mutations
GB1_doped_mut = unique(GB1_doped_mut)
#get rid of mutations that don't change codons and those that give stop codons
GB1_doped_mut = GB1_doped_mut[WTaa != Mut & Mut != "*"]
print(paste("unique possible NT1  mutations, w/o PTC:",nrow(GB1_doped_mut)))
print(paste0("fraction unique possible NNS mutations (",55*19,"): ",round(nrow(GB1_doped_mut)/55/19*100,digits=1),"%"))
GB1_doped_mut[,PosMut := paste0(Pos,Mut)]

#create doped double data.tables
#full dataset
doubles = fread("GB1/processed_data/DMS_doubles_preE.txt")
doubles[,':=' (PosMut1 = paste0(Pos1,Mut1),PosMut2 = paste0(Pos2,Mut2))]
doubles = doubles[PosMut1 %in% GB1_doped_mut$PosMut & PosMut2 %in% GB1_doped_mut$PosMut]
write.table(x = doubles, file = paste0(dataset_dir,"processed_data/DMS_doubles_preE_doped.txt"),
            quote = F,row.names = F, col.names = T)

#from read-downsampled versions
for (RD in read_downsampling) {
  doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles_preE_RD",RD,".txt"))
  doubles[,':=' (PosMut1 = paste0(Pos1,Mut1),PosMut2 = paste0(Pos2,Mut2))]
  doubles = doubles[PosMut1 %in% GB1_doped_mut$PosMut & PosMut2 %in% GB1_doped_mut$PosMut]
  write.table(x = doubles, file = paste0(dataset_dir,"processed_data/DMS_doubles_preE_doped_RD",RD,".txt"),
              quote = F,row.names = F, col.names = T)
}

##################################################################
### calculate epistasis null model and call pos./neg.epistasis ###
##################################################################

#read all downsampled data files
double_files = list.files(path = paste0(dataset_dir,"processed_data/"))[grep("DMS_doubles_preE_",list.files(path = paste0(dataset_dir,"processed_data/")))]
ID = sapply(double_files,FUN=function(X){gsub(".txt","",gsub("DMS_doubles_preE_","",X))})

for (idx in seq_along(ID)) {
  doubles = fread(paste0(dataset_dir,"processed_data/",double_files[idx]))
  doubles2 = copy(doubles[,.SD,,.SDcols = c(1:16)])
  #lower bound of fitness
  lower_bound_F = doubles2[is.fitness == T & is.reads0 == T & fitness1 + fitness2 < -8,median(fitness,na.rm=T)]
  #call epistatic interactions
  doubles = call_epistasis_binary(doubles,
                                  lower_bound_F,
                                  dataset_dir = dataset_dir,
                                  output_filename = paste0("DMS_doubles_",ID[idx],".txt"),
                                  prefix = paste0("GB1_",ID[idx]))
}


#############################################
### calculate pairwise interaction scores ###
#############################################

double_files = list.files(path = paste0(dataset_dir,"processed_data/"))[grep("DMS_doubles_preE_",list.files(path = paste0(dataset_dir,"processed_data/")))]
ID = sapply(double_files,FUN=function(X){gsub(".txt","",gsub("DMS_doubles_preE_","",X))})
for (idx in seq_along(ID)) {
  doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles_",ID[idx],".txt"))
  PWI = calculate_pairwise_interaction_scores(doubles,
                                              dataset_dir = dataset_dir,
                                              output_filename = paste0("DMS_PWI_",ID[idx],".txt"),
                                              detailed = F)
}

#### assemble all the combined scores from all PWI data.tables into one
PWI_complete = fread("GB1/processed_data/DMS_PWI.txt")[,.(Pos1,Pos2,WT_AA1,WT_AA2,combined_score)]
for (idx in seq_along(ID)) {
  PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_",ID[idx],".txt"))
  PWI_complete = merge(PWI_complete,
                       PWI[,.(Pos1,Pos2,V1=combined_score)],
                       by=c("Pos1","Pos2"),all=T)
  names(PWI_complete)[ncol(PWI_complete)] = paste0("combined_score_",ID[idx])
}
write.table(PWI_complete,file=paste0(dataset_dir,"processed_data/DMS_PWI_complete.txt"),
            row.names=F,quote = F)


#########################################
### deep contact transform PWI scores ###
#########################################

prefix = "GB1_"

PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_complete.txt"))
PWI_transformed = deepcontact_transform_basic2d(PWI[,-grep("WT_AA",names(PWI)),with=F],
                                                dataset_dir = dataset_dir,
                                                output_filename = "DMS_PWI_complete_deepcontact.txt",
                                                deepcontact_dir = "where/is/deepcontact/",
                                                prefix = prefix)




###########################################################
######## predict secondary structure from PWI data ########
###########################################################

PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_complete.txt"))
#### predict secondary structure
predict_secondary_structure_elements(PWI,
                                     dataset_dir = dataset_dir,
                                     prefix = prefix,
                                     known_SS = "GB1/processed_data/PDB_secondary_structure_1pga_A.txt")

#### predict beta sheets
predict_beta_sheets(PWI,
                    input_ss0 = fread(paste0(dataset_dir,"processed_data/",prefix,"_secondary_structure_prediction.txt")),
                    dataset_dir = dataset_dir,
                    prefix = prefix,
                    known_ss_file = "GB1/processed_data/PDB_secondary_structure_1pga_A.txt",
                    known_bsi_file = "GB1/processed_data/PDB_beta_sheet_hbonds_1pga.txt")




##########################################################################################
##### evaluate predicted contacts (top scoring pairs) against reference structure ########
##########################################################################################

contactmap = fread("GB1/processed_data/PDB_contactmap_1pga_A.txt")

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


#########################################
##### XPLOR structure prediction ########
#########################################

#use predicted tertiary contacts and secondary structure elements to generate structural models

#### use  DMS data to derive restraints for tertiary contacts & use secondary structure restraints derived from PSIPRED predictions (no beta sheet pairing restraints)
# +compare to deepcontact transformed scores

PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI_complete.txt"))
PWI_DC = fread(paste0(dataset_dir,"processed_data/DMS_PWI_complete_deepcontact.txt"))
names(PWI_DC)[names(PWI_DC) %in% setdiff(names(PWI_DC),c("Pos1","Pos2","WT_AA1","WT_AA2"))] = paste0(setdiff(names(PWI_DC),c("Pos1","Pos2","WT_AA1","WT_AA2")),"_DC")
PWI2 = merge(PWI,PWI_DC,by=c("Pos1","Pos2"),all=T)

XPLOR_wrapper(input_PWI = PWI2,
              SS_mode = "SSonly",
              input_SS_file = "PSIPRED_secondary_structure.txt", ### copy this from GB1 folder into the dataset folder
              prefix = "GB1_ds_PSIPRED_SSonly_",
              dataset_dir = dataset_dir,
              cores = 15,
              queue = "short-sl7,long-sl7",
              protein_sequence = scan("GB1/dataset/GB1_sequence.fasta",what = "character")[2],
              pdb_file = "GB1/dataset/PDB/g_xray.pdb",
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

analyse_XPLOR_results(XPLOR_dir = paste0(dataset_dir,"results/XPLOR/GB1_ds_PSIPRED_SSonly_/"),
                      contactmap = fread("GB1/processed_data/contactmap_1pga_A.txt"),
                      draw_contactmaps = F)
