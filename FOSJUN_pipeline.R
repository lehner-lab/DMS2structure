##############################
##### GB1 from Olson2014 #####
##############################

#this is the pipeline used to analyse FOS/JUN data from Diss&Lehner 2018

#first, set the working directory to the DMS2structure folder
setwd("/where/is/DMS2structure/")


#source scripts
filelist = list.files('scripts/')
sapply(paste0('scripts/',filelist),source,.GlobalEnv)

#create the necessary subfolder structure for all results and processed data
dataset_dir = "FOSJUN/"
create_directory_structure(dataset_dir)
#then save this script in the dataset_dir
#and paste all necessary source data into dataset_dir/dataset

#load required packages
require(data.table)
require(ggplot2)
require(cowplot)
require(GGally)


#############################################################################################
##### preprocess data (calculate fitness scores and errors, set quality thresholds etc) #####
#############################################################################################

##### first download the raw count table from here
# https://www.dropbox.com/s/7h5swkpxlk6sleh/Diss2018_trans_Q20.txt?dl=0
# and save it to the following path:
output = fread(paste0(dataset_dir,"dataset/Diss2018_trans_Q20.txt"))
output[wt_aa_F == mut_aa_F,':=' (aa_pos_F = NA,mut_aa_F = NA)]
output[wt_aa_J == mut_aa_J,':=' (mut_aa_J = NA,aa_pos_J = NA)]

output2 = output[,.(count_r1_t0 = sum(i1),count_r2_t0 = sum(i2),count_r3_t0 = sum(i3),
          count_r1_t1 = sum(o1),count_r2_t1 = sum(o2),count_r3_t1 = sum(o3)),
       .(aa_pos_F,aa_pos_J,mut_aa_F,mut_aa_J)]
names(output2)[1:4] = c("Pos1","Pos2","Mut1","Mut2")
output2[,Pos1 := Pos1 + 1] #correct positions starting from 0
output2[,Pos2 := Pos2 + 1]
output2[,FOSid := paste0(Pos1,Mut1),.(Pos1,Mut1)]
output2[,JUNid := paste0(Pos2,Mut2),.(Pos2,Mut2)]
output2[,Nmut := as.numeric(!is.na(Mut1)) + as.numeric(!is.na(Mut2)),.(FOSid,JUNid)]

## plot replicate scores for singles
ggpairs(output2[Nmut==1][sample(.N,1000),log10(.SD+1),.SDcols = grep("count",names(output2))])
ggsave(file = paste0(dataset_dir,"results/preprocessing/all_vars_Nmut1.pdf"),width=9,height=9)
# and doubles
ggpairs(output2[Nmut==2][sample(.N,1000),log10(.SD+1),.SDcols = grep("count",names(output2))])
ggsave(file = paste0(dataset_dir,"results/preprocessing/all_vars_Nmut2.pdf"),width=9,height=9)

##### get rid of lower read count input peak in double mutants
all_data = output2[count_r1_t0>=10 & count_r2_t0>=10 & count_r3_t0>=10]

#replot
ggpairs(all_data[Nmut==2][sample(.N,1000),log10(.SD+1),.SDcols = grep("count",names(all_data))])
ggsave(file = paste0(dataset_dir,"results/preprocessing/all_vars_10readcountthreshold_Nmut2.pdf"),width=9,height=9)

##### mutant statistics
output2[,.N,Nmut]
all_data[,.N,Nmut]

##### fitness scores for all three replicates from counts for all variants
wildtype_helper = all_data[Nmut==0]
all_data[,fitness_r1 := log(count_r1_t1/count_r1_t0/(wildtype_helper$count_r1_t1/wildtype_helper$count_r1_t0))]
all_data[,fitness_r2 := log(count_r2_t1/count_r2_t0/(wildtype_helper$count_r2_t1/wildtype_helper$count_r2_t0))]
all_data[,fitness_r3 := log(count_r3_t1/count_r3_t0/(wildtype_helper$count_r3_t1/wildtype_helper$count_r3_t0))]

##### Poissonian errors from counts
all_data[,sigma_r1 := sqrt(1/count_r1_t1 + 1/count_r1_t0 + 1/wildtype_helper$count_r1_t1 + 1/wildtype_helper$count_r1_t0)]
all_data[,sigma_r2 := sqrt(1/count_r2_t1 + 1/count_r2_t0 + 1/wildtype_helper$count_r2_t1 + 1/wildtype_helper$count_r2_t0)]
all_data[,sigma_r3 := sqrt(1/count_r3_t1 + 1/count_r3_t0 + 1/wildtype_helper$count_r3_t1 + 1/wildtype_helper$count_r3_t0)]

###### build wildtype,singles and doubles data.tables
wildtype = all_data[Nmut==0]
singles = rbind(all_data[Nmut==1 & !is.na(Pos1),cbind(Pos = Pos1, Mut = Mut1,protein = "FOS",id = paste0("FOS_",Pos1,Mut1),
                                                      .SD[,1:6],count_t0_mean = rowMeans(.SD[,7:9])),,
                         .SDcols = c(grep("fitness",names(all_data)),grep("sigma",names(all_data)),grep("count_r[123]_t0",names(all_data)))],
                all_data[Nmut==1 & !is.na(Pos2),cbind(Pos = Pos2, Mut = Mut2,protein = "JUN",id = paste0("JUN_",Pos2,Mut2),
                                                      .SD[,1:6],count_t0_mean = rowMeans(.SD[,7:9])),,
                         .SDcols = c(grep("fitness",names(all_data)),grep("sigma",names(all_data)),grep("count_r[123]_t0",names(all_data)))])
doubles = all_data[Nmut==2,cbind(Pos1,Pos2,Mut1,Mut2,id1 = paste0("FOS_",Pos1,Mut1),id2 = paste0("JUN_",Pos2,Mut2),.SD),,
                   .SDcols = c(grep("count",names(all_data)),grep("fitness",names(all_data)),grep("sigma",names(all_data)))]




########################################################################
##### bayesian framework for fitness estimation for double mutants #####
########################################################################

doubles[,bin_count_r1_t0 := findInterval(log10(count_r1_t0),seq(0.5,4,0.25))]
doubles[,.(.N,mean(count_r1_t0)),bin_count_r1_t0][order(bin_count_r1_t0)]
ggplot(doubles[bin_count_r1_t0 < 10],aes(fitness_r1,..count..,color=factor(bin_count_r1_t0))) +
  geom_density(adjust=1)
#>> low coverage leads to systematic bias at low fitness values
#>> estimate what fitness scores are for variants with low sequence coverage

#precomputed data (computation takes ~30min)
doubles = fread(paste0(dataset_dir,"processed_data/doubles_bayesian_cond.txt"))

############################################################################################################
#1) poisson distribution for score likelihood
lam_d = 0.025

## calculate posterior double mutant fitness based on prior from single mutants
postpois_conditioned_singleF = function(i){
  require(data.table)

  count_in = double_data[i,count_t0]
  count_out = double_data[i,count_t1]
  lam_in = exp(seq(floor(log(count_in+0.1)-max(c(0.5,1/log10(count_in+1.75)))),(log(count_in+0.1)+max(c(0.5,1/log10(count_in+1.75)))),lam_d))
  lam_out = exp(seq(floor(log(count_out+0.1)-max(c(0.5,1/log10(count_out+1.75)))),(log(count_out+0.1)+max(c(0.5,1/log10(count_out+1.75)))),lam_d))
  lam_low = range(log(lam_out))[1] - range(log(lam_in))[2]
  lam_high = range(log(lam_out))[2] - range(log(lam_in))[1]
  idx = row(matrix(NA,nrow=length(lam_out),ncol=length(lam_in))) - col(matrix(NA,nrow=length(lam_out),ncol=length(lam_in)))
  likelihood = sapply(split(outer(dpois(count_out,lambda = lam_out),dpois(count_in,lambda = lam_in)),idx),sum)
  score_prior = density(score_prior_cond[,.(fdist = sqrt((double_data[i,fitness1]-fitness1)^2+(double_data[i,fitness2]-fitness2)^2),fitness)][
                                  order(fdist)][1:Nneighbours,fitness],
                        from = (lam_low-wt_corr),
                        to = (lam_high-wt_corr),
                        n = as.integer(as.character(round((lam_high-lam_low)/lam_d + 1)))) #super weird bug
  posterior = score_prior$y*likelihood

  moments = list()
  moments[1] = weighted.mean(x = score_prior$x,w = posterior)
  moments[2] = sqrt(sum(( moments[[1]]-score_prior$x)^2 * posterior)/
                      sum(posterior))
  return(moments)
}


require(parallel)
# Use the detectCores() function to find the number of cores in system
no_cores <- detectCores() - 1
# Setup cluster
clust <- makeCluster(no_cores) #This line will take time

for (rep in 1:3) {
  #wildtype "correction" to calculate scores
  wt_corr = unlist(wildtype[,log(.SD[,2]/.SD[,1]),.SDcols = grep(paste0("count_r",rep),names(wildtype))])
  Nneighbours = 1000
  #data for prior calculation
  double_idx = unlist(doubles[,unlist(.SD[,1]) > 100 & unlist(.SD[,2]) > -Inf,,.SDcols = c(paste0("count_r",rep,"_t0"),paste0("fitness_r",rep))])
  double_data = doubles[,.(id1,id2,count_t0 = get(paste0("count_r",rep,"_t0")),count_t1 = get(paste0("count_r",rep,"_t1")),
                           fitness = get(paste0("fitness_r",rep)))]
  double_data[,fitness1 := singles[id == id1,get(paste0("fitness_r",rep))],id1]
  double_data[,fitness2 := singles[id == id2,get(paste0("fitness_r",rep))],id2]
  score_prior_cond = double_data[double_idx]

  # make variables available to each core's workspace
  clusterExport(clust, list("double_data","lam_d","wt_corr","score_prior_cond","Nneighbours"))

  #posterior fitness conditioned on single fitness
  t=proc.time()
  helper = parSapply(clust,X = 1:nrow(double_data), postpois_conditioned_singleF)
  print(proc.time()-t)
  helper1 = matrix(unlist(helper),nrow=2)
  doubles[,paste0("fitness_cond_r",rep) := helper1[1,]]
  doubles[,paste0("sigma_cond_r",rep) := helper1[2,]]
}

stopCluster(clust)
# write estimates to file for later reload in case of rerunning the script
write.table(x = doubles,file = paste0(dataset_dir,"processed_data/doubles_bayesian_cond.txt"),
            quote = F, row.names = F, col.names = T)

##### plot results from fitness estimations
helper = melt(doubles[count_r1_t0 >= 1,.SD,,.SDcols = c(grep("^count_r1_t0",names(doubles)),grep("fitness.*_r1",names(doubles)))],id.vars = "count_r1_t0")
p1 = ggplot(helper,aes(count_r1_t0,value)) +
  geom_hex() +
  scale_fill_distiller(direction = 1,trans="log10") +
  scale_x_log10() +
  facet_grid(. ~ variable)

require(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

helper = melt(doubles[between(bin_count_r1_t0,1,9),.SD,,
                      .SDcols = c(grep("bin_count",names(doubles)),grep("fitness.*_r1",names(doubles)))],id.vars = "bin_count_r1_t0")
numcolors = length(unique(helper$bin_count_r1_t0))
p2 = ggplot(helper,aes(value,..density..,color = factor(bin_count_r1_t0))) +
  geom_density() +
  scale_color_manual(values=getPalette(numcolors)) +
  facet_grid(. ~ variable) +
  labs(color="bin")

helper = rbind(doubles[count_r1_t0 >= 1,.(fitness= fitness_r1,sigma = sigma_r1,type = "uncorr")],
               doubles[count_r1_t0 >= 1,.(fitness= fitness_cond_r1,sigma = sigma_cond_r1,type = "cond")])
helper[,type := factor(type,levels = c("uncorr","cond"))]
p3 = ggplot(helper,aes(fitness,sigma)) +
  geom_hex() +
  scale_y_log10() +
  scale_fill_distiller(direction = 1,trans="log10") +
  facet_grid(. ~ type)

P = plot_grid(plotlist = list(p1,p2,p3),nrow=3)
P
ggsave(plot = P,file = paste0(dataset_dir,"results/preprocessing/fitness_calculation_uncorrVScond.pdf"),width=6,height=8)
############################################################################################################

##############################################################
##### rescale fitness values with Guillaume's parameters #####
##############################################################
##### relation between ppi scores and "fitness" scores
# Guillaume normalized the fitness over the whole experiment by the number of generations yeast cells have grown
#1) get ppi scores from supplementary table 1
output3 = fread(paste0(dataset_dir,"dataset/elife-32472-supp1-v2.txt"))
all_data = merge(all_data,output3[,.(Pos1=pos1,Pos2=pos2,Mut1=mut1,Mut2=mut2,
                                          ppi_r1 = do_ppi1, ppi_r2 = do_ppi2, ppi_r3 = do_ppi3)],
                 by=c("Pos1","Pos2","Mut1","Mut2"),all.x=T)

ggplot(all_data,aes(fitness_r1,ppi_r1)) +
  geom_point() +
  geom_smooth(method="lm")
#linear relationship between score and ppi --> fit for conversion
conv_factors = list()
conv_factors[[1]] = all_data[,lm(ppi_r1 ~ fitness_r1)]
conv_factors[[2]] = all_data[,lm(ppi_r2 ~ fitness_r2)]
conv_factors[[3]] = all_data[,lm(ppi_r3 ~ fitness_r3)]

##### transform uncorrected fitness values
fitness_idx = names(singles)[grep("fitness",names(singles))]
sigma_idx = names(singles)[grep("sigma",names(singles))]
for (idx in seq_along(fitness_idx)) {
  #log fitness scores 
  singles[,fitness_idx[idx] := log(.SD * conv_factors[[idx]]$coefficients[2] + conv_factors[[idx]]$coefficients[1]),.SDcols = fitness_idx[idx]]
  doubles[,fitness_idx[idx] := log(.SD * conv_factors[[idx]]$coefficients[2] + conv_factors[[idx]]$coefficients[1]),.SDcols = fitness_idx[idx]]
  #relative errors
  wildtype[,sigma_idx[idx] := .SD[,1] * conv_factors[[idx]]$coefficients[2] / exp(.SD[,2]),,.SDcols = c(sigma_idx[idx],fitness_idx[idx])]
  singles[,sigma_idx[idx] := .SD[,1] * conv_factors[[idx]]$coefficients[2] / exp(.SD[,2]),,.SDcols = c(sigma_idx[idx],fitness_idx[idx])]
  doubles[,sigma_idx[idx] := .SD[,1] * conv_factors[[idx]]$coefficients[2] / exp(.SD[,2]),,.SDcols = c(sigma_idx[idx],fitness_idx[idx])]
} #ignore warnings
##### transform conditioned fitness values
fitness_idx = names(doubles)[grep("fitness_cond",names(doubles))]
sigma_idx = names(doubles)[grep("sigma_cond",names(doubles))]
for (idx in seq_along(fitness_idx)) {
  #log fitness scores 
  doubles[,fitness_idx[idx] := log(.SD * conv_factors[[idx]]$coefficients[2] + conv_factors[[idx]]$coefficients[1]),.SDcols = fitness_idx[idx]]
  #relative errors
  doubles[,sigma_idx[idx] := .SD[,1] * conv_factors[[idx]]$coefficients[2] / exp(.SD[,2]),,.SDcols = c(sigma_idx[idx],fitness_idx[idx])]
}

############################
##### merge replicates #####
############################

idx_vec = c("","_cond")
for (idx in 1:2) {
  f_vec = doubles[,grep(paste0("fitness",idx_vec[idx]),names(doubles)),with=F]
  s_vec = doubles[,grep(paste0("sigma",idx_vec[idx]),names(doubles)),with=F]
  
  doubles[,paste0("fitness",idx_vec[idx]) := rowSums(f_vec/(s_vec^2),na.rm=T) / rowSums(1/(s_vec^2),na.rm=T)]
  doubles[,paste0("sigma",idx_vec[idx]) := sqrt(1 / rowSums(1/(s_vec^2),na.rm=T))]
}

ggpairs(doubles[sample(.N,1000),cbind(.SD[,1:2],log10(.SD[,3:4])),,
                .SDcols = c("fitness","fitness_cond","sigma","sigma_cond")],
        aes(alpha=0.1))
ggsave(paste0(dataset_dir,"results/preprocessing/fitness_merged_uncorrVScond.pdf"),width=7.5,height=7.5)

### also merge wildtype and singles
wildtype[,fitness := 0]
wildtype[,sigma := sqrt(1/rowSums(1/.SD^2)),,.SDcols = grep("sigma",names(wildtype))]

singles[,fitness := rowSums(.SD[,1:3]/(.SD[,4:6]^2),na.rm=T)/
          rowSums(1/(.SD[,4:6]^2),na.rm=T),,
        .SDcols = c(grep("fitness",names(singles)),grep("sigma",names(singles)))]
singles[,sigma := sqrt(1/rowSums(1/(.SD[,1:3]^2),na.rm=T)),,
        .SDcols = grep("sigma",names(singles))]

########################
## finalize data.tables
singles = singles[!is.na(fitness) & !is.na(sigma)]

doubles[,count_t0_mean := rowMeans(.SD),,.SDcols = c("count_r1_t0","count_r2_t0","count_r3_t0")]
doubles = doubles[,setdiff(1:ncol(doubles),grep("r[123]",names(doubles))),with=F]

doubles[,fitness1 := singles[id == unique(unlist(.SD)),fitness],by=.(id1),.SDcols = "id1"]
doubles[,sigma1 := singles[id == unique(unlist(.SD)),sigma],by=.(id1),.SDcols = "id1"]
doubles[,fitness2 := singles[id == unique(unlist(.SD)),fitness],by=.(id2),.SDcols = "id2"]
doubles[,sigma2 := singles[id == unique(unlist(.SD)),sigma],by=.(id2),.SDcols = "id2"]
doubles = doubles[!is.na(fitness1) & !is.na(fitness2) & !is.na(sigma1) & !is.na(sigma2)]

###### mark STOP variants
singles[Mut=="*",STOP := protein]
doubles[Mut1 == "*" ,STOP:="FOS"]
doubles[Mut2 == "*" ,STOP:="JUN"]
doubles[Mut1 == "*" & Mut2 == "*" ,STOP:="FOSJUN"]

#and compare them
helper = rbind(singles[,.(fitness,count_t0_mean,STOP,mut="single",type = "uncorr")],
               doubles[,.(fitness,count_t0_mean,STOP,mut="double",type = "uncorr")],
               singles[,.(fitness,count_t0_mean,STOP,mut="single",type = "cond")],
               doubles[,.(fitness=fitness_cond,count_t0_mean,STOP,mut="double",type = "cond")])
helper[,type := factor(type,levels=c("uncorr","cond"))]
p1 = ggplot(helper[!is.na(STOP)],aes(y=count_t0_mean,x=fitness,color=STOP,shape = mut,..scaled..)) +
  geom_density2d() +
  scale_y_log10() +
  coord_cartesian(xlim=c(-1.2,0.1)) +
  facet_grid(. ~ type)
p2 = ggplot(helper,aes(fitness,color=STOP,linetype = mut,..scaled..)) + 
  geom_density(adjust=0.75) +
  coord_cartesian(xlim=c(-1.2,0.1)) +
  facet_grid(. ~ type)
plot_grid(plotlist = list(p1,p2),nrow=2)
# singles with a STOP mutation have higher fitness than doubles with TWO STOP mutations
# doubles with only one STOP mutants are bimodal
ggsave(paste0(dataset_dir,"results/preprocessing/fitness_density_STOPvars.pdf"),width=8,height=7)




###### mark variants with nonsensical fitness values
doubles[,is.fitness:=F]
doubles[is.na(STOP),is.fitness:=T]

###### define which variants have enough read coverage
doubles[,is.reads0:=T]
# this should be (count_t0_mean > 100) for uncorr fitness

###### add wild-type amino acid to tables
WT_F = unique(output[,.(aa_pos_F,wt_aa_F)])[!is.na(aa_pos_F)][order(aa_pos_F),.(Pos0 = aa_pos_F+1,WT_AA = wt_aa_F)]
WT_J = unique(output[,.(aa_pos_J,wt_aa_J)])[!is.na(aa_pos_J)][order(aa_pos_J),.(Pos0 = aa_pos_J+1,WT_AA = wt_aa_J)]
singles[protein=="FOS",WT_AA := WT_F[Pos0 == Pos,WT_AA],.(Pos)]
singles[protein=="JUN",WT_AA := WT_J[Pos0 == Pos,WT_AA],.(Pos)]
doubles[,WT_AA1 := WT_F[Pos0 == Pos1,WT_AA],Pos1]
doubles[,WT_AA2 := WT_J[Pos0 == Pos2,WT_AA],Pos2]

###### estimate lower measurement limit of fitness assay 
xd1=density((singles[,fitness]),bw=.02)
plot(xd1)
xd2=density((doubles[is.fitness==T & !is.na(fitness_cond),fitness_cond]),bw=.02)
plot(xd2)
lb = xd1$x[xd1$y==max(xd1$y[xd1$x < -0.5])] #-.75

ggplot(doubles,aes(fitness1+fitness2,fitness_cond)) +
  geom_hex() +
  scale_fill_continuous(trans="log10") +
  geom_abline() +
  geom_smooth() +
  geom_hline(yintercept = lb) +
  geom_hline(yintercept = doubles[fitness1 + fitness2 < -1.2 & is.fitness==T,weighted.mean(fitness_cond,w=1/sigma_cond^2)],color="red") +
  geom_hline(yintercept = doubles[STOP=="FOS",weighted.mean(fitness_cond,w=1/sigma_cond^2)],color="green",linetype=2) +
  geom_hline(yintercept = doubles[STOP=="JUN",weighted.mean(fitness_cond,w=1/sigma_cond^2)],color="green",linetype=2) +
  geom_hline(yintercept = doubles[STOP=="FOSJUN",weighted.mean(fitness_cond,w=1/sigma_cond^2)],color="darkgreen",linetype=2)
ggsave(paste0(dataset_dir,"results/preprocessing/fitness_expVSobs_lowerfitnesslimit.pdf"),width=8,height=7)

#lower bound estimate from double dead variants and double STOPs coindicde very well, use this as lower bond
lower_bound_F = doubles[fitness1 + fitness2 < -1.2 & is.fitness==T,weighted.mean(fitness_cond,w=1/sigma_cond^2)]
print(lower_bound_F)


###### reorder doubles table
doubles = doubles[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,id1,id2,STOP,is.fitness,is.reads0,count_t0_mean,
                     fitness1,sigma1,fitness2,sigma2,
                     fitness,sigma,
                     fitness_cond,sigma_cond)]


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
#only for conditioned fitness estimates
doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles_preE.txt"))
doubles2 = copy(doubles[,.SD,,.SDcols = c(1:16,grep(names(doubles),pattern="_cond$"))])
names(doubles2)[c(17,18)] = c("fitness","sigma")
doubles = call_epistasis_binary(doubles2,
                                lower_bound_F,
                                dataset_dir = dataset_dir,
                                prefix = "FOSJUN_",
                                sym=F)

#############################################
### calculate pairwise interaction scores ###
#############################################

# doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles.txt"))
PWI = calculate_pairwise_interaction_scores(doubles,
                                            N_resample = 10^2,
                                            dataset_dir = dataset_dir,
                                            modus = "trans",
                                            detailed = F)


#########################################
### deep contact transform PWI scores ###
#########################################

# PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))
PWI_transformed = deepcontact_transform_basic2d(PWI[,.(Pos1,Pos2,
                                                       epistasis_score)],
                                                dataset_dir = dataset_dir,
                                                deepcontact_dir = "where/is/deepcontact/",
                                                modus = "trans",
                                                prefix = "FOSJUN_")

### negative controls for DeepContact
# 3x permutate epistasis_score, while keeping matrix symmetry
set.seed(1603)
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))[,.(Pos1,Pos2,control1 = sample(epistasis_score),
                                                                 control2 = sample(epistasis_score),
                                                                 control3 = sample(epistasis_score))]
write.table(x = PWI, file = paste0(dataset_dir,"processed_data/DMS_PWI_DC_control.txt"),
            quote = F,row.names = F, col.names = T)
PWI_transformed = deepcontact_transform_basic2d(PWI,
                                                dataset_dir = dataset_dir,
                                                output_filename = "DMS_PWI_DC_control_deepcontact.txt",
                                                deepcontact_dir = "where/is/deepcontact/",
                                                modus = "trans",
                                                prefix = "FOSJUN_")


##################################################################################
##### read distance and secodnary structure info from PDB files (or PSIPRED) #####
##################################################################################

#deposit PDB files in the dataset subfolders first

singles = fread(paste0(dataset_dir,"processed_data/DMS_singles.txt"))
WT_seq_chainE = paste0(unique(singles[grep(singles$id,pattern="FOS"),.(Pos,WT_AA)])[order(Pos)][,WT_AA],collapse="")
WT_seq_chainF = paste0(unique(singles[grep(singles$id,pattern="JUN"),.(Pos,WT_AA)])[order(Pos)][,WT_AA],collapse="")

#calculate contact maps for trans interaction
pairdistances_from_PDB(input_file = paste0(dataset_dir,"dataset/PDB/1fos_1.pdb"),
                       dataset_dir = dataset_dir,
                       given_chainids = c("E","F"),
                       aa_seq = list(WT_seq_chainE,WT_seq_chainF),
                       idx_pdb_start=c(162,286),
                       idx_DMS_start = c(1,1),
                       debug_this=F)

#for FOS cis interactions
pairdistances_from_PDB(input_file = paste0(dataset_dir,"dataset/PDB/1fos_1.pdb"),
                       dataset_dir = dataset_dir,
                       given_chainids = "E",
                       aa_seq = WT_seq_chainE,
                       idx_pdb_start=162,
                       idx_DMS_start = 1,
                       debug_this=F)

#for JUN cis interactions
pairdistances_from_PDB(input_file = paste0(dataset_dir,"dataset/PDB/1fos_1.pdb"),
                       dataset_dir = dataset_dir,
                       given_chainids = "F",
                       aa_seq = WT_seq_chainE,
                       idx_pdb_start= 286,
                       idx_DMS_start = 1,
                       debug_this=F)


###########################################
###### evaluate epistasis data (QC) #######
###########################################
### these are some basic scripts to evaluate the dataset and predictions

doubles = fread(paste0(dataset_dir,"processed_data/DMS_doubles.txt"))
PWI = fread(paste0(dataset_dir,"processed_data/DMS_PWI.txt"))
prefix = "FOSJUN_"
contactmap = fread(paste0(dataset_dir,"processed_data/PDB_contactmap_1fos_1_EF.txt"))
contactmap_FOS = fread(paste0(dataset_dir,"processed_data/PDB_contactmap_1fos_1_E.txt"))
contactmap_JUN = fread(paste0(dataset_dir,"processed_data/PDB_contactmap_1fos_1_F.txt"))

#check how subsets over which positive or negative epistasis can be classified are distributed in single mutant fitness space
epistasis_analytics_subsets_singlemutantspace(doubles,
                                              modus = "trans",
                                              dataset_dir = dataset_dir,prefix = prefix)

#marginal distribution of # of variants suitable for epistasis classification across position pairs
epistasis_analytics_NumEvars_marginal(doubles,
                                      modus = "trans",
                                      dataset_dir = dataset_dir,prefix = prefix)

#variants suitable for epistasis classification versus median fitness of single mutants
epistasis_analytics_NumEvars_fitness(doubles,
                                     modus = "trans",
                                     dataset_dir = dataset_dir,prefix = prefix)

#check the spatial distribution of the number of variants suitable for epistasis classification
epistasis_analytics_NumEvars_spatial(PWI,
                                     modus = "trans",
                                     dataset_dir = dataset_dir,prefix = prefix)

# cumulative distribution function of epistatic variants as function of distance 
# (if a known contactmap is available)
epistasis_analytics_subsets_CDF(doubles,
                                dataset_dir = dataset_dir,
                                contactmap = contactmap,
                                prefix = prefix,
                                modus = "trans",
                                dist_type = "scHAmin", lindist = 5)



###########################################################
######## predict secondary structure from PWI data ########
###########################################################

#for FOS
predict_secondary_structure_elements(PWI[,.(Pos1,Pos2,
                                            association_score1)],
                                     dataset_dir = dataset_dir,
                                     prefix = "FOS_",
                                     known_SS = paste0(dataset_dir,"processed_data/PDB_secondary_structure_1fos_1_E.txt"),
                                     scale_long = 1/6^2)

#for JUN
predict_secondary_structure_elements(PWI[,.(Pos1,Pos2,
                                            association_score2)],
                                     dataset_dir = dataset_dir,
                                     prefix = "JUN_",
                                     known_SS = paste0(dataset_dir,"processed_data/PDB_secondary_structure_1fos_1_F.txt"),
                                     scale_long = 1/6^2)


##########################################################################################
##### evaluate predicted contacts (top scoring pairs) against reference structure ########
##########################################################################################

#### true positive rate of top contacts + contactmaps + eCDFs
evaluate_contacts_vs_PDB(contacts = PWI[,.(Pos1,Pos2,epistasis_score)],
                         contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                         secondary_structure=NA,
                         dataset_dir = dataset_dir,
                         lindist=5,
                         modus = "trans",
                         prefix = prefix)

#### interaction scores versus distance in reference structure
score_vs_distance_scatter(contacts = PWI[,.(Pos1,Pos2,epistasis_score)],
                          contactmap = contactmap[,.(Pos1,Pos2,scHAmin)],
                          dataset_dir = dataset_dir,
                          modus = "trans",
                          prefix = prefix)

