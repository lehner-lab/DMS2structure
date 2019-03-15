
############################################################
### calculate pairwise interaction scores via resampling ###
############################################################

calculate_pairwise_interaction_scores = function(
  double_data,
  N_resample = 10^4,
  modus = "cis",
  dataset_dir,
  output_filename = "DMS_PWI.txt",
  diagonal_entries = "one",
  detailed = F,
  cores = NULL) {
  
  ### variables 
  # double_data: the doubles data.table
  # N_resample: number of resampling runs
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # output_filename: filename to write datatable to in dataset_dir/processed_data/
  # modus: "cis" (single protein) or "trans" (protein-protein interaction)
  # diagonal_entries: what the diagonal entries should be when filling up symmetric enrichment matrices
  #                   - one: sets diagonal to 1 (default); however, in data with very sparse epistatic interactions, this will give a high bias for poisition pairs with high enrichments to also have high correlations
  #                   - means: use mean enrichment over all pairs a position is involved in
  # detailed: if FALSE, it will only give epistasis, association and combined scores as outputs; 
  #           if TRUE, it will output all the intermediate scores and their uncertainities; see below in the script for naming conventions
  # cores: number of cores to request on the computing cluster

  require(data.table)
  require(stringr)
  require(gdata)
  require(corpcor)
  
  #data needed from doubles data.table
  A = copy(double_data[is.fitness==TRUE&is.reads0==TRUE,
                       .(Pos1,Pos2,WT_AA1,WT_AA2,fitness,sigmaE,pos_epistasis,neg_epistasis,F_fit_upper,F_fit_lower)])
  
  #positions in protein
  pos_range = c(min(c(A$Pos1,A$Pos2)),max(c(A$Pos1,A$Pos2)))
  #initialize pairwise interaction datatable
  PWI = data.table(Pos1=rep(pos_range[1]:pos_range[2],times=diff(pos_range)+1),
                   Pos2=rep(pos_range[1]:pos_range[2],each=diff(pos_range)+1))
  #add wild-type amino acids and number of variants available for positive or negative epistasis classification
  PWI = merge(PWI,
              unique(A[,.(WT_AA1,WT_AA2,NposE = sum(pos_epistasis),NnegE = sum(neg_epistasis)),.(Pos1,Pos2)]),
              by=c("Pos1","Pos2"),all.x=T)
  setkey(PWI,Pos1,Pos2)
  
  ## calculate enrichment priors
  pseq = seq(0,1,0.001)
  # calculate prior of positive epistasis enrichment as density of measured enrichments
  prior_posE = density(x=A[pos_epistasis==TRUE,
                           .(.N,posE = sum(pnorm(F_fit_upper,fitness,sigmaE,lower.tail = F),na.rm=T)/.N),
                           by=.(Pos1,Pos2)][N>max(c(median(N),10)),posE],from=0, to=1, n=length(pseq),bw=0.05)
  
  # calculate prior of negative epistasis enrichment as density of measured enrichments
  prior_negE = density(x=A[neg_epistasis==TRUE,
                           .(.N,negE = sum(pnorm(F_fit_lower,fitness,sigmaE,lower.tail = T),na.rm=T)/.N),
                           by=.(Pos1,Pos2)][N>max(c(median(N),10)),negE],from=0, to=1, n=length(pseq),bw=0.05)
  #plot 
  ggplot(data=rbind(data.frame(x=prior_posE$x,y=prior_posE$y,type="prior_posE"),
                    data.frame(x=prior_negE$x,y=prior_negE$y,type="prior_negE")),
         aes(x,y,color=type)) + 
    geom_line() +
    labs(x="fraction enriched",y="density")
  ggsave(paste0(dataset_dir,"results/preprocessing/enrichment_prior_posE_negE.pdf"),height=3.5,width=4.5)
  
  ### calculate/resample interaction scores in parallel
  require(parallel)
  no_cores <- detectCores()-1
  if(!is.null(cores)){
    no_cores <- cores
  }
  clust <- makeCluster(no_cores)
  
  ######### positive epistasis enrichment
  posE_enr = function(i) {
    require(data.table)
    require(gdata)
    ### draw fitness values according to error
    DT = A[,fhelper := rnorm(.N,fitness,sigmaE)][pos_epistasis==T,.(NposE = sum(fhelper>F_fit_upper),Ntotal = .N),.(Pos1,Pos2)]
    setkey(DT,Pos1,Pos2)
    DT[,posE :=  pseq[which.min(abs(cumsum(dbinom(NposE, size=Ntotal, prob=pseq)*prior_posE$y)/
                                      sum(dbinom(NposE, size=Ntotal, prob=pseq)*prior_posE$y)-runif(1)))],.(Pos1,Pos2)]
    ## transfer enrichments to matrix
    posE_enr_matrix = matrix(NA,nrow=diff(pos_range)+1,ncol=diff(pos_range)+1)
    posE_enr_matrix[cbind(DT$Pos1-pos_range[1]+1,DT$Pos2-pos_range[1]+1)] = DT$posE
    ## fill na values
    if (modus == "cis") {
      #only upper triangle has valid data
      posE_enr_matrix[is.na(posE_enr_matrix) & upper.tri(posE_enr_matrix)] = sample(c(posE_enr_matrix[!is.na(posE_enr_matrix)]),sum(is.na(posE_enr_matrix) & upper.tri(posE_enr_matrix)),replace = T)  
      #make symmetric
      lowerTriangle(posE_enr_matrix) = upperTriangle(posE_enr_matrix,byrow=T)  
    } else { #modus == "trans"
      #full matrix has valid data
      posE_enr_matrix[is.na(posE_enr_matrix)] = sample(c(posE_enr_matrix[!is.na(posE_enr_matrix)]),sum(is.na(posE_enr_matrix)),replace = T)
    }
    #output linearized matrix
    return(c(posE_enr_matrix))
  }
  clusterSetRNGStream(cl = clust,1234567)
  clusterExport(clust, list("A","prior_posE","pseq","modus","pos_range"),envir = environment())
  posE_enr_rs = parSapply(clust,1:N_resample,posE_enr) #this is a position X #resampling matrix
  
  
  ######### correlation of positive epistasis patterns
  posE_cor = function(i) {
    require(corpcor)
    ## add "pseudocount" for logging
    posE_enr_matrix = matrix(posE_enr_rs[,i],nrow=diff(pos_range)+1) + quantile(posE_enr_rs[,i],probs = 0.25,na.rm=T)
    
    if (modus == "cis") {
      if (diagonal_entries == "one") {## set diagonal to 1
        diag(posE_enr_matrix) = 1
      } else if (diagonal_entries ==  "means") {
        diag(posE_enr_matrix) = rowMeans(posE_enr_matrix,na.rm=T)
      }
    }
    ### calculate correlation on log values
    if (do.transpose == T) {
      return(c(corpcor::cor.shrink(t(log(posE_enr_matrix)),verbose=F)))  
    } else {
      return(c(corpcor::cor.shrink(log(posE_enr_matrix),verbose=F)))
    }
  }
  if (modus == "cis") { #if these are cis-interactions with a symmetric data matrix, dimension over which correlations are calculate doesn't matter
    do.transpose = F
    posE_cor_rs = sapply(1:N_resample,posE_cor)
  } else { #if this is a trans-interaction dataset with a non-symmetric data matrix, 
    # the dimensions over which the correlations are calculate reveal different info about the underlying proteins
    #### NOTE: currently only works if both proteins (/molecules) have the same amout of positions mutated
    #protein 1
    do.transpose = T
    posE1_cor_rs = sapply(1:N_resample,posE_cor)
    #protein 2
    do.transpose = F
    posE2_cor_rs = sapply(1:N_resample,posE_cor)
  }
  
  ######### partial correlation of positive epistasis patterns
  posE_pcor = function(i) {
    require(corpcor)
    ## add "pseudocount" for logging
    posE_enr_matrix = matrix(posE_enr_rs[,i],nrow=diff(pos_range)+1) + quantile(posE_enr_rs[,i],probs = 0.25,na.rm=T)
    if (modus == "cis") {
      if (diagonal_entries == "one") {## set diagonal to 1
        diag(posE_enr_matrix) = 1
      } else if (diagonal_entries ==  "means") {
        diag(posE_enr_matrix) = rowMeans(posE_enr_matrix,na.rm=T)
      }
    }
    ### calculate correlation on log values
    if (do.transpose == T) {
      return(c(corpcor::pcor.shrink(t(log(posE_enr_matrix)),verbose=F)))  
    } else {
      return(c(corpcor::pcor.shrink(log(posE_enr_matrix),verbose=F)))
    }
  }
  if (modus == "cis") {
    do.transpose = F
    posE_pcor_rs = sapply(1:N_resample,posE_pcor)
  } else {
    do.transpose = T
    posE1_pcor_rs = sapply(1:N_resample,posE_pcor)
    
    do.transpose = F
    posE2_pcor_rs = sapply(1:N_resample,posE_pcor)
  }
  
  
  ######### negative epistasis enrichment
  negE_enr = function(i) {
    require(data.table)
    require(gdata)
    ### draw fitness values according to error
    DT = A[,fhelper := rnorm(.N,fitness,sigmaE)][neg_epistasis==T,.(NnegE = sum(fhelper<F_fit_lower),Ntotal = .N),.(Pos1,Pos2)]
    setkey(DT,Pos1,Pos2)
    DT[,negE :=  pseq[which.min(abs(cumsum(dbinom(NnegE, size=Ntotal, prob=pseq)*prior_negE$y)/
                                      sum(dbinom(NnegE, size=Ntotal, prob=pseq)*prior_negE$y)-runif(1)))],.(Pos1,Pos2)]
    ## transfer enrichments to matrix
    negE_enr_matrix = matrix(NA,nrow=diff(pos_range)+1,ncol=diff(pos_range)+1)
    negE_enr_matrix[cbind(DT$Pos1-pos_range[1]+1,DT$Pos2-pos_range[1]+1)] = DT$negE
    ## fill na values
    if (modus == "cis") {
      negE_enr_matrix[is.na(negE_enr_matrix) & upper.tri(negE_enr_matrix)] = sample(c(negE_enr_matrix[!is.na(negE_enr_matrix)]),sum(is.na(negE_enr_matrix) & upper.tri(negE_enr_matrix)),replace = T)
      lowerTriangle(negE_enr_matrix) = upperTriangle(negE_enr_matrix,byrow=T)  
    } else {
      negE_enr_matrix[is.na(negE_enr_matrix)] = sample(c(negE_enr_matrix[!is.na(negE_enr_matrix)]),sum(is.na(negE_enr_matrix)),replace = T)
    }
    return(c(negE_enr_matrix))
  }
  clusterSetRNGStream(cl = clust,1234567)
  clusterExport(clust, list("A","prior_negE","pseq","modus","pos_range"),envir = environment())
  negE_enr_rs = parSapply(clust,1:N_resample,negE_enr)
  
  ######### correlation of negative epistasis patterns
  negE_cor = function(i) {
    require(corpcor)
    ## add "pseudocount" for logging
    negE_enr_matrix = matrix(negE_enr_rs[,i],nrow=diff(pos_range)+1) + quantile(negE_enr_rs[,i],probs = 0.25,na.rm=T)
    if (modus == "cis") {# set diagonal to 1
      if (diagonal_entries == "one") {## set diagonal to 1
        diag(negE_enr_matrix) = 1
      } else if (diagonal_entries ==  "means") {
        diag(negE_enr_matrix) = rowMeans(negE_enr_matrix,na.rm=T)
      }
    }
    ### calculate correlation on log values
    if (do.transpose == T) {
      return(c(corpcor::cor.shrink(t(log(negE_enr_matrix)),verbose=F)))  
    } else {
      return(c(corpcor::cor.shrink(log(negE_enr_matrix),verbose=F)))
    }
  }
  if (modus == "cis") {
    do.transpose = F
    negE_cor_rs = sapply(1:N_resample,negE_cor)
  } else {
    #protein 1
    do.transpose = T
    negE1_cor_rs = sapply(1:N_resample,negE_cor)
    #protein 2
    do.transpose = F
    negE2_cor_rs = sapply(1:N_resample,negE_cor)
  }
  
  ######### partial correlation of negative epistasis patterns
  negE_pcor = function(i) {
    require(corpcor)
    ## add "pseudocount" for logging
    negE_enr_matrix = matrix(negE_enr_rs[,i],nrow=diff(pos_range)+1) + quantile(negE_enr_rs[,i],probs = 0.25,na.rm=T)
    if (modus == "cis") {
      if (diagonal_entries == "one") {## set diagonal to 1
        diag(negE_enr_matrix) = 1
      } else if (diagonal_entries ==  "means") {
        diag(negE_enr_matrix) = rowMeans(negE_enr_matrix,na.rm=T)
      }
    }
    ### calculate correlation on log values
    if (do.transpose == T) {
      return(c(corpcor::pcor.shrink(t(log(negE_enr_matrix)),verbose=F)))  
    } else {
      return(c(corpcor::pcor.shrink(log(negE_enr_matrix),verbose=F)))
    }
  }
  if (modus == "cis") {
    do.transpose = F
    negE_pcor_rs = sapply(1:N_resample,negE_pcor)
  } else {
    #protein 1
    do.transpose = T
    negE1_pcor_rs = sapply(1:N_resample,negE_pcor)
    #protein 2
    do.transpose = F
    negE2_pcor_rs = sapply(1:N_resample,negE_pcor)
  }
  
  stopCluster(clust)
  
  ### transfer enrichments and partial correlations to data.table
  setkey(PWI,Pos1,Pos2)
  #enrichments
  #posE
  PWI[,posE_enr_mean := apply(posE_enr_rs,1,mean)]
  PWI[,posE_enr_sd := apply(posE_enr_rs,1,sd)]
  PWI[,posE_enr := posE_enr_mean/posE_enr_sd]
  
  #negE
  PWI[,negE_enr_mean := apply(negE_enr_rs,1,mean)]
  PWI[,negE_enr_sd := apply(negE_enr_rs,1,sd)]
  PWI[,negE_enr := negE_enr_mean/negE_enr_sd]
  
  
  #(partial) correlations
  if (modus == "cis") {
    PWI[,posE_cor_mean := apply(posE_cor_rs,1,mean)]
    PWI[,posE_cor_sd := apply(posE_cor_rs,1,sd)]
    PWI[,posE_cor := posE_cor_mean/apply(posE_cor_rs,1,sd)]
    PWI[,posE_pcor_mean := apply(posE_pcor_rs,1,mean)]
    PWI[,posE_pcor_sd := apply(posE_pcor_rs,1,sd)]
    PWI[,posE_pcor := posE_pcor_mean/apply(posE_pcor_rs,1,sd)]
    PWI[,negE_cor_mean := apply(negE_cor_rs,1,mean)]
    PWI[,negE_cor_sd := apply(negE_cor_rs,1,sd)]
    PWI[,negE_cor := negE_cor_mean/apply(negE_cor_rs,1,sd)]
    PWI[,negE_pcor_mean := apply(negE_pcor_rs,1,mean)]
    PWI[,negE_pcor_sd := apply(negE_pcor_rs,1,sd)]
    PWI[,negE_pcor := negE_pcor_mean/apply(negE_pcor_rs,1,sd)]
  } else { #modus == "trans"
    PWI[,posE1_cor_mean := apply(posE1_cor_rs,1,mean)]
    PWI[,posE1_cor_sd := apply(posE1_cor_rs,1,sd)]
    PWI[,posE1_cor := posE1_cor_mean/apply(posE1_cor_rs,1,sd)]
    PWI[,posE2_cor_mean := apply(posE2_cor_rs,1,mean)]
    PWI[,posE2_cor_sd := apply(posE2_cor_rs,1,sd)]
    PWI[,posE2_cor := posE2_cor_mean/apply(posE2_cor_rs,1,sd)]
    
    PWI[,posE1_pcor_mean := apply(posE1_pcor_rs,1,mean)]
    PWI[,posE1_pcor_sd := apply(posE1_pcor_rs,1,sd)]
    PWI[,posE1_pcor := posE1_pcor_mean/apply(posE1_pcor_rs,1,sd)]
    PWI[,posE2_pcor_mean := apply(posE2_pcor_rs,1,mean)]
    PWI[,posE2_pcor_sd := apply(posE2_pcor_rs,1,sd)]
    PWI[,posE2_pcor := posE2_pcor_mean/apply(posE2_pcor_rs,1,sd)]
    
    PWI[,negE1_cor_mean := apply(negE1_cor_rs,1,mean)]
    PWI[,negE1_cor_sd := apply(negE1_cor_rs,1,sd)]
    PWI[,negE1_cor := negE1_cor_mean/apply(negE1_cor_rs,1,sd)]
    PWI[,negE2_cor_mean := apply(negE2_cor_rs,1,mean)]
    PWI[,negE2_cor_sd := apply(negE2_cor_rs,1,sd)]
    PWI[,negE2_cor := negE2_cor_mean/apply(negE2_cor_rs,1,sd)]
    
    PWI[,negE1_pcor_mean := apply(negE1_pcor_rs,1,mean)]
    PWI[,negE1_pcor_sd := apply(negE1_pcor_rs,1,sd)]
    PWI[,negE1_pcor := negE1_pcor_mean/apply(negE1_pcor_rs,1,sd)]
    PWI[,negE2_pcor_mean := apply(negE2_pcor_rs,1,mean)]
    PWI[,negE2_pcor_sd := apply(negE2_pcor_rs,1,sd)]
    PWI[,negE2_pcor := negE2_pcor_mean/apply(negE2_pcor_rs,1,sd)]
  }
  
  #### calculate epistasis score  by weighted averaging over positive and negative enrichments
  PWI[,posE_negE_enr_mean := sum(c((posE_enr_mean/(posE_enr_sd^2 )),(negE_enr_mean/(negE_enr_sd^2 ))),na.rm=T)/
        sum(c((1/posE_enr_sd^2 ),(1/negE_enr_sd^2 )),na.rm=T),.(Pos1,Pos2)]
  PWI[,posE_negE_enr_sd := sqrt(1/sum(c((1/posE_enr_sd^2),(1/negE_enr_sd^2)),na.rm=T)),.(Pos1,Pos2)]
  PWI[,epistasis_score := posE_negE_enr_mean/posE_negE_enr_sd]
  
  #### calculate association scores (and merged correlation score) by weighted averaging over positive and negative (parital) correlations
  if (modus == "cis") {
    PWI[,posE_negE_pcor_mean := (posE_pcor_mean/apply(posE_pcor_rs,1,sd)^2 + negE_pcor_mean/apply(negE_pcor_rs,1,sd)^2)/(1/apply(posE_pcor_rs,1,sd)^2 + 1/apply(negE_pcor_rs,1,sd)^2)]
    PWI[,posE_negE_pcor_sd := sqrt(1/(1/apply(posE_pcor_rs,1,sd)^2 + 1/apply(negE_pcor_rs,1,sd)^2))]
    PWI[,association_score := posE_negE_pcor_mean/sqrt(1/(1/apply(posE_pcor_rs,1,sd)^2 + 1/apply(negE_pcor_rs,1,sd)^2))]
    PWI[,posE_negE_cor_mean := (posE_cor_mean/apply(posE_cor_rs,1,sd)^2 + negE_cor_mean/apply(negE_cor_rs,1,sd)^2)/(1/apply(posE_cor_rs,1,sd)^2 + 1/apply(negE_cor_rs,1,sd)^2)]
    PWI[,posE_negE_cor_sd := sqrt(1/(1/apply(posE_cor_rs,1,sd)^2 + 1/apply(negE_cor_rs,1,sd)^2))]
    PWI[,posE_negE_cor := posE_negE_cor_mean/sqrt(1/(1/apply(posE_cor_rs,1,sd)^2 + 1/apply(negE_cor_rs,1,sd)^2))] 
  } else { #modus == "trans"
    PWI[,posE1_negE1_pcor_mean := (posE1_pcor_mean/apply(posE1_pcor_rs,1,sd)^2 + negE1_pcor_mean/apply(negE1_pcor_rs,1,sd)^2)/(1/apply(posE1_pcor_rs,1,sd)^2 + 1/apply(negE1_pcor_rs,1,sd)^2)]
    PWI[,posE1_negE1_pcor_sd := sqrt(1/(1/apply(posE1_pcor_rs,1,sd)^2 + 1/apply(negE1_pcor_rs,1,sd)^2))]
    PWI[,association_score1 := posE1_negE1_pcor_mean/sqrt(1/(1/apply(posE1_pcor_rs,1,sd)^2 + 1/apply(negE1_pcor_rs,1,sd)^2))]
    PWI[,posE2_negE2_pcor_mean := (posE2_pcor_mean/apply(posE2_pcor_rs,1,sd)^2 + negE2_pcor_mean/apply(negE2_pcor_rs,1,sd)^2)/(1/apply(posE2_pcor_rs,1,sd)^2 + 1/apply(negE2_pcor_rs,1,sd)^2)]
    PWI[,posE2_negE2_pcor_sd := sqrt(1/(1/apply(posE2_pcor_rs,1,sd)^2 + 1/apply(negE2_pcor_rs,1,sd)^2))]
    PWI[,posE2_negE2_pcor := posE2_negE2_pcor_mean/sqrt(1/(1/apply(posE2_pcor_rs,1,sd)^2 + 1/apply(negE2_pcor_rs,1,sd)^2))]
    
    PWI[,posE1_negE1_cor_mean := (posE1_cor_mean/apply(posE1_cor_rs,1,sd)^2 + negE1_cor_mean/apply(negE1_cor_rs,1,sd)^2)/(1/apply(posE1_cor_rs,1,sd)^2 + 1/apply(negE1_cor_rs,1,sd)^2)]
    PWI[,posE1_negE1_cor_sd := sqrt(1/(1/apply(posE1_cor_rs,1,sd)^2 + 1/apply(negE1_cor_rs,1,sd)^2))]
    PWI[,association_score2 := posE1_negE1_cor_mean/sqrt(1/(1/apply(posE1_cor_rs,1,sd)^2 + 1/apply(negE1_cor_rs,1,sd)^2))]
    PWI[,posE2_negE2_cor_mean := (posE2_cor_mean/apply(posE2_cor_rs,1,sd)^2 + negE2_cor_mean/apply(negE2_cor_rs,1,sd)^2)/(1/apply(posE2_cor_rs,1,sd)^2 + 1/apply(negE2_cor_rs,1,sd)^2)]
    PWI[,posE2_negE2_cor_sd := sqrt(1/(1/apply(posE2_cor_rs,1,sd)^2 + 1/apply(negE2_cor_rs,1,sd)^2))]
    PWI[,posE2_negE2_cor := posE2_negE2_cor_mean/sqrt(1/(1/apply(posE2_cor_rs,1,sd)^2 + 1/apply(negE2_cor_rs,1,sd)^2))]
    
  }
  
  ### calculate combined score
  if (modus == "cis") { #only works for cis-interations within a protein
    
    #using all epistasis info
    PWI[Pos1 != Pos2,temp1 := scale(association_score)]
    PWI[Pos1 != Pos2,temp2 := scale(epistasis_score)]
    PWI[Pos1 != Pos2,combined_score := sum(c(temp1,temp2),na.rm=T),.(Pos1,Pos2)]
    
    #only using positive epistasis info
    PWI[,c("temp1","temp2") := NULL]
    PWI[Pos1 != Pos2,temp1 := scale(posE_pcor)]
    PWI[Pos1 != Pos2,temp2 := scale(posE_enr)]
    PWI[Pos1 != Pos2,combined_score_posE := sum(c(temp1,temp2),na.rm=T),.(Pos1,Pos2)]
    PWI[,c("temp1","temp2") := NULL]
    
    #only using negative epistasis info
    PWI[Pos1 != Pos2,temp1 := scale(negE_pcor)]
    PWI[Pos1 != Pos2,temp2 := scale(negE_enr)]
    PWI[Pos1 != Pos2,combined_score_negE := sum(c(temp1,temp2),na.rm=T),.(Pos1,Pos2)]
    PWI[,c("temp1","temp2") := NULL]
  }
  
  if (detailed) {
    PWI2 = PWI
  } else {
    if (modus=="cis") {
      PWI2 = PWI[,.(Pos1,Pos2,WT_AA1,WT_AA2,NposE,NnegE,
                    epistasis_score,
                    association_score,
                    combined_score)]
    } else {
      PWI2 = PWI[,.(Pos1,Pos2,WT_AA1,WT_AA2,NposE,NnegE,
                    epistasis_score,
                    association_score1,
                    association_score2)]
    }
  }
  
  ## write data.table to file
  write.table(x = PWI2, file = paste0(dataset_dir,"processed_data/",output_filename),
              quote = F,row.names = F, col.names = T)
  
  return(PWI2)
}



