######################################################################
####### plot epistasis subset regions in single fitness space ########
######################################################################
epistasis_analytics_subsets_singlemutantspace = function(doubles,
                                                         dataset_dir,
                                                         prefix = "",
                                                         modus = "cis") {
  
  ### variables 
  # doubles: doubles data.table
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # modus: "cis" (single protein) or "trans" (protein-protein interaction)
  
  if (modus == "trans") {
    data = copy(doubles[is.fitness==TRUE & is.reads0==TRUE])
  } else { #if cis library make symmetric
    data = switch_double_DT(doubles[is.fitness==TRUE & is.reads0==TRUE],list(c("fitness1","fitness2")),c("fitness","pos_epistasis","neg_epistasis"))  
  }
  
  ## bin data in both fitness directions
  x1 = seq(data[,quantile(fitness1,0,na.rm=T)*0.99],data[,quantile(fitness1,1,na.rm=T)*1.01],length.out = 25)
  data[,bin1:=findInterval(fitness1,x1)]
  data[,bin2:=findInterval(fitness2,x1)]
  data[,X1 := x1[bin1],by=bin1]
  data[,X2 := x1[bin2],by=bin2]
  
  ## 2d histogram
  setkey(data,X1,X2)
  N_binned_all = data[.(rep(x1,length(x1)),(rep(x1,each=length(x1)))),.(N_all=.N),by=.EACHI]
  N_binned_pos = data[pos_epistasis==TRUE][.(rep(x1,length(x1)),(rep(x1,each=length(x1)))),.(N_pos=.N),by=.EACHI]
  N_binned_neg = data[neg_epistasis==TRUE][.(rep(x1,length(x1)),(rep(x1,each=length(x1)))),.(N_neg=.N),by=.EACHI]
  
  ## smooth minimally to facilitate contour plot drawing
  for (i in 1:nrow(N_binned_all)) {
    w = exp(-((N_binned_all[,X1]-N_binned_all[i,X1])^2 + (N_binned_all[,X2]-N_binned_all[i,X2])^2)/0.001) 
    set(N_binned_all,i,"N_all_smooth",sum(N_binned_all$N_all*w)/sum(w))
    set(N_binned_pos,i,"N_pos_smooth",sum(N_binned_pos$N_pos*w)/sum(w))
    set(N_binned_neg,i,"N_neg_smooth",sum(N_binned_neg$N_neg*w)/sum(w))
  }
  
  N_binned_merged = merge(merge(N_binned_all,N_binned_pos,by=c("X1","X2"),all = TRUE),N_binned_neg,by=c("X1","X2"),all = TRUE)
  N_binned_merged[N_all_smooth < 1,N_all_smooth := 0]
  N_binned_merged[,N_pos_smooth := as.numeric(N_pos_smooth >= 1)]
  N_binned_merged[,N_neg_smooth := as.numeric(N_neg_smooth >= 1)]
  
  ggplot() +
    geom_raster(data=N_binned_merged,aes(X1,X2,fill=N_all),interpolate=F) +
    scale_fill_gradient(low="gray75",high="dodgerblue4",trans="log10",na.value = "white") +
    geom_contour(data=N_binned_merged,aes(X1,X2,z=N_neg_smooth,color="gold"),size=1,bins=1,alpha=0.75) +
    geom_contour(data=N_binned_merged,aes(X1,X2,z=N_pos_smooth,color="red"),size=1,bins=1,alpha=0.75) +
    scale_color_manual(values=c("gold","red"),labels=c("negative","positive")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x="single mutant fitness 1",y="single mutant fitness 2",
         color="epistasis subsets",fill="#variants")
  ggsave(paste0(dataset_dir,"results/epistasis/",prefix,"epistasis_subsets_singlemutantspace.pdf"),height=4,width=6)
}


##########################################################################
####### plot marginal distribution of the number of variants suitable ####
####### for epistasis classification over all position pairs #############
##########################################################################
epistasis_analytics_NumEvars_marginal = function(doubles,
                                                 dataset_dir,
                                                 prefix = "",
                                                 modus = "cis") {
  ### variables 
  # doubles: doubles data.table
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # modus: "cis" (single protein) or "trans" (protein-protein interaction)
  
  theme_set(theme_classic())
  
  if (modus == "cis") {
    doubles_sym = switch_double_DT(doubles[is.reads0==T & is.fitness == T],list(c("Pos1","Pos2"),c("fitness1","fitness2")),c("pos_epistasis","neg_epistasis"))  
  } else {
    doubles_sym = copy(doubles[is.reads0==T & is.fitness == T])
  }
  
  DT_numbervars = doubles_sym[,.(num_all = .N,
                                 num_posE = sum(pos_epistasis==T),
                                 num_negE = sum(neg_epistasis==T)),
                              .(Pos1,Pos2)]
  
  DT_numbervars_melt = melt(DT_numbervars,id.vars = "Pos1",measure.vars = c("num_all","num_posE","num_negE"))

  ggplot(DT_numbervars_melt,aes(value,color=variable,..count..)) +
    geom_density(adjust=0.5) +
    scale_color_manual(breaks=c("num_all","num_posE","num_negE"),values = c("black","red","gold"),
                       labels = c(paste0("all, <n> = ",DT_numbervars[,round(mean(num_all))]),
                                  paste0("pos.E, <n> = ",DT_numbervars[,round(mean(num_posE))]),
                                  paste0("neg.E, <n> = ",DT_numbervars[,round(mean(num_negE))]))) +
    scale_x_continuous(limits = c(0,361),breaks = seq(0,350,50),expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0,350,50),expand = c(0,0)) +
    labs(x="number of double mutants per position pair",y="density [a.u.]",
         color = "data subset")
  ggsave(paste0(dataset_dir,"results/epistasis/",prefix,"number_epistatic_variants.pdf"),width=5,height=3)
}

####################################################################
####### number of variants suitable for epistasis classification ###
#######  versus single mutant fitness ##############################
epistasis_analytics_NumEvars_fitness = function(doubles,
                                                dataset_dir,
                                                prefix = "",
                                                modus = "cis") {
  
  ### variables 
  # doubles: doubles data.table
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # modus: "cis" (single protein) or "trans" (protein-protein interaction)
  
  theme_set(theme_classic(base_size = 9))
  
  if (modus == "cis") {
    doubles_sym = switch_double_DT(doubles[is.reads0==T & is.fitness == T],list(c("Pos1","Pos2"),c("fitness1","fitness2")),c("pos_epistasis","neg_epistasis"))  
    DT_numbervars = doubles_sym[,.(num_all = .N,
                                   num_posE = sum(pos_epistasis==T),
                                   num_negE = sum(neg_epistasis==T),
                                   median_fitness = median(fitness1,na.rm=T),
                                   variants = "all"),
                                .(Pos1,Pos2)]
    ##### fitness versus number of variants
    DT_mN_mF = DT_numbervars[,.(median_N = median(num_all),median_fitness = mean(median_fitness),variants = "all"),Pos1]
    DT_mN_mF = rbind(DT_mN_mF,DT_numbervars[,.(median_N = median(num_posE),median_fitness = mean(median_fitness),variants = "posE"),Pos1])
    DT_mN_mF = rbind(DT_mN_mF,DT_numbervars[,.(median_N = median(num_negE),median_fitness = mean(median_fitness),variants = "negE"),Pos1])
    DT_mN_mF[,variants := factor(variants,levels = c("all","posE","negE"))]
    
    ggplot(DT_mN_mF,aes(median_fitness,median_N,color = variants)) +
      geom_point() +
      scale_color_manual(values = c("black","red","gold")) +
      geom_smooth(se=F) +
      labs(x = "median single mutant fitness at position",
           y = "median # of double mutants in pairs involving position")
    ggsave(paste0(dataset_dir,"results/epistasis/",prefix,"fitness_vs_numberdoublemutants.pdf"),width=5,height=4)
    
  } else {
    # doubles_sym = copy(doubles[is.reads0==T & is.fitness == T])
    DT_numbervars = doubles[is.reads0==T & is.fitness == T,
                            .(num_all = .N,
                              num_posE = sum(pos_epistasis==T),
                              num_negE = sum(neg_epistasis==T),
                              median_fitness1 = median(fitness1,na.rm=T),
                              median_fitness2 = median(fitness2,na.rm=T),
                              variants = "all"),
                            .(Pos1,Pos2)]
    ##### fitness versus number of variants
    DT_mN_mF1 = DT_numbervars[,.(median_N = median(num_all),median_fitness = mean(median_fitness1),variants = "all"),Pos1]
    DT_mN_mF1 = rbind(DT_mN_mF1,DT_numbervars[,.(median_N = median(num_posE),median_fitness = mean(median_fitness1),variants = "posE"),Pos1])
    DT_mN_mF1 = rbind(DT_mN_mF1,DT_numbervars[,.(median_N = median(num_negE),median_fitness = mean(median_fitness1),variants = "negE"),Pos1])
    DT_mN_mF1[,variants := factor(variants,levels = c("all","posE","negE"))]
    
    ggplot(DT_mN_mF1,aes(median_fitness,median_N,color = variants)) +
      geom_point() +
      scale_color_manual(values = c("black","red","gold")) +
      geom_smooth(se=F) +
      labs(x = "median single mutant fitness at position",
           y = "median number of double mutants in position pairs involving position")
    ggsave(paste0(dataset_dir,"results/epistasis/",prefix,"fitness_vs_numberdoublemutants_protein1.pdf"),width=5,height=4)
    
    DT_mN_mF2 = DT_numbervars[,.(median_N = median(num_all),median_fitness = mean(median_fitness2),variants = "all"),Pos2]
    DT_mN_mF2 = rbind(DT_mN_mF2,DT_numbervars[,.(median_N = median(num_posE),median_fitness = mean(median_fitness2),variants = "posE"),Pos2])
    DT_mN_mF2 = rbind(DT_mN_mF2,DT_numbervars[,.(median_N = median(num_negE),median_fitness = mean(median_fitness2),variants = "negE"),Pos2])
    DT_mN_mF2[,variants := factor(variants,levels = c("all","posE","negE"))]
    ggplot(DT_mN_mF2,aes(median_fitness,median_N,color = variants)) +
      geom_point() +
      scale_color_manual(values = c("black","red","gold")) +
      geom_smooth(se=F) +
      labs(x = "median single mutant fitness at position",
           y = "median number of double mutants in position pairs involving position")
    ggsave(paste0(dataset_dir,"results/epistasis/",prefix,"fitness_vs_numberdoublemutants_protein2.pdf"),width=5,height=4)
  }
}


###############################################################
##### spatial distribution of variants per position pair ######
###############################################################
epistasis_analytics_NumEvars_spatial = function(PWI,
                                                dataset_dir,
                                                prefix = "",
                                                modus = "cis") {
  
  ### variables 
  # PWI: pairwise interaction score data.table
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # modus: "cis" (single protein) or "trans" (protein-protein interaction)
  
  if (modus == "cis") { #plot positive and negative subsets as halves of same matrix
    ggplot() +
      geom_raster(data = PWI[Pos1<Pos2],aes(Pos1,Pos2,fill=NposE)) +
      geom_raster(data = PWI[Pos1<Pos2],aes(Pos2,Pos1,fill=NnegE)) +
      scale_fill_distiller(direction=1) +
      scale_x_continuous(expand = c(0.01,0),breaks = seq(5,55,5)) +
      scale_y_reverse(expand = c(0.01,0),breaks = seq(5,55,5)) +
      labs(x = "position",y="position",title = "posE: lower left, negE: upper right",fill = "# variants")
    ggsave(paste0(dataset_dir,"results/epistasis/",prefix,"number_variants_per_positionpair.pdf"),width=5.5,height=4.2)
  } else { #separate matrices
    Ppos = ggplot() +
      geom_raster(data = PWI,aes(Pos1,Pos2,fill=NposE)) +
      scale_fill_distiller(direction=1,limits = c(0,250)) +
      scale_x_continuous(expand = c(0.01,0),breaks = seq(5,55,5)) +
      scale_y_reverse(expand = c(0.01,0),breaks = seq(5,55,5)) +
      labs(x = "protein 1",y="protein 2",title = "posE variants",fill = "# variants")
    
    Pneg = ggplot() +
      geom_raster(data = PWI,aes(Pos2,Pos1,fill=NnegE)) +
      scale_fill_distiller(direction=1,limits = c(0,250)) +
      scale_x_continuous(expand = c(0.01,0),breaks = seq(5,55,5)) +
      scale_y_reverse(expand = c(0.01,0),breaks = seq(5,55,5)) +
      labs(x = "protein 1",y="protein 2",title = "negE variants",fill = "# variants")
    
    plot_grid(plotlist = list(Ppos,Pneg),nrow=1)
    
    ggsave(paste0(dataset_dir,"results/epistasis/",prefix,"number_variants_per_positionpair.pdf"),width=8,height=3.2)
    
  }
}






######################################################################
####### CDF of epistasis variants versus distance (if contactmap avail.)
######################################################################
epistasis_analytics_subsets_CDF = function(doubles,
                                           dataset_dir,
                                           prefix = "",
                                           contactmap,
                                           modus = "cis",
                                           dist_type = "scHAmin",
                                           dist_cutoff = 8,
                                           lindist = 5) { 
  
  ### variables
  # doubles: doubles data.table
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # contactmap: known pair-wise distance contactmap from a PDB file
  # modus: "cis" (single protein) or "trans" (protein-protein interaction)
  # dist_type: which distance metric to use
  # dist_cutoff: if not NA, will be used to report the fraction of variants below that distances
  # lindist: only use position pairs with a sequence separation greater than this (only for cis libraries)
  
  data = merge(doubles[is.fitness==TRUE & is.reads0==TRUE],contactmap,by=c("Pos1","Pos2"))

  if (modus == "trans") {
    data1 = rbind(data[,c(.SD,type="all"),,.SDcols = dist_type], #all variants that pass read filter threshold
                  data[pos_epistasis==T,c(.SD,type="posEsubset"),,.SDcols = dist_type], #all variants suitable for positive epistasis classification
                  data[pos_epistasis==T & fitness > F_fit_upper,c(.SD,type="posE"),,.SDcols = dist_type], #in this subset, variants that are classified as positive epistatic
                  data[neg_epistasis==T,c(.SD,type="negEsubset"),,.SDcols = dist_type], #same for negative
                  data[neg_epistasis==T & fitness < F_fit_lower,c(.SD,type="negE"),,.SDcols = dist_type])
    title_label = "trans"
  } else {
    data1 = rbind(data[Pos1<Pos2-lindist,c(.SD,type="all"),,.SDcols = dist_type], #all variants that pass read filter threshold
                  data[Pos1<Pos2-lindist & pos_epistasis==T,c(.SD,type="posEsubset"),,.SDcols = dist_type],#all variants suitable for positive epistasis classification
                  data[Pos1<Pos2-lindist & pos_epistasis==T & fitness > F_fit_upper,c(.SD,type="posE"),,.SDcols = dist_type],#in this subset, variants that are classified as positive epistatic
                  data[Pos1<Pos2-lindist & neg_epistasis==T,c(.SD,type="negEsubset"),,.SDcols = dist_type],#same for negative
                  data[Pos1<Pos2-lindist & neg_epistasis==T & fitness < F_fit_lower,c(.SD,type="negE"),,.SDcols = dist_type])
    title_label = paste0("cis, lindist > ",lindist)
  }
  names(data1)[1] = "distance"
  
  data1[,type := factor(type,levels=c("all","posE","posEsubset","negE","negEsubset"))]
  setkey(data1,type)
  ggplot(data1,aes(x=distance,color=type)) + 
    stat_ecdf() +
    scale_color_manual(values = c("black","red","orange","gold","yellow"),
                       labels = data1[,.(paste0(type," ",round(sum(.SD<dist_cutoff)/.N,digits=2)*100,"%")),type,.SDcols = "distance"]$V1) +
    scale_x_continuous(expand=c(0,0),breaks = seq(0,50,5)) +
    scale_y_continuous(expand=c(0,0),breaks = seq(0,1,0.25)) +
    geom_vline(xintercept = dist_cutoff,linetype=2) +
    labs(x=paste0(dist_type," distance [A]"),y="cumulative probability",title = title_label) +
    theme_classic()
  ggsave(ifelse(modus == "trans",
                paste0(dataset_dir,"results/epistasis/",prefix,"epistasis_subsets_CDF_",dist_type,".pdf"),
                paste0(dataset_dir,"results/epistasis/",prefix,"epistasis_subsets_CDF_",dist_type,"_lindist",lindist,".pdf")),height=3,width=4)
}


