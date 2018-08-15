#####################################################################
### calculate non-parametric epistasis null model, call epistasis & #
#####  define subsets for positive/neative epistasis evaluation #####
#####################################################################

call_epistasis_binary = function(double_data,
                                 lower_bound_F,
                                 dataset_dir,
                                 output_filename = "DMS_doubles.txt",
                                 prefix = "",
                                 xsig=2,
                                 sym = T,
                                 Q = 0.05,
                                 epistasis_error_from_slopes = T) {
  
  ### variables 
  # double_data: the doubles data.table
  # lower_bound_F: an estimate of the lower bound of the fitness assay
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/epistasis/
  # output_filename: filename to write datatable to in dataset_dir/processed_data/
  # prefix:to be added to results files (in case of running diff. versions of data from same dataset etc)
  # xsig: significance threshold for calling epistasis, 2 is fine
  # sym: logical for whether the position-position fitness map is symmetrical (only FALSE for protein-protein interactions!!)
  # Q: the percentile used for upper (1-Q) and lower (Q) fitness surface calculation
  # epistasis_error_from_slopes: logical, should sigmeE (error for epistasis values) be calculate by taking slopes of median fitness surface into account
  
  set.seed(1603)
  
  require(data.table)
  require(mgcv)
  require(caTools)
  
  ##################################################################################
  ######## calculate quantile fitness surfaces for epistasis classification ########
  ##################################################################################
  
  #make position-position fitness map symmetrical
  if (sym == T) {
    DT = switch_double_DT(double_data[is.fitness == TRUE & is.reads0 == TRUE],
                          cols_switchdouble = list(c("fitness1","fitness2"),c("sigma1","sigma2")),
                          cols_double = c("fitness","sigma"))
  } else { #for non-symmetrical protein-protein interaction leave as is
    DT = copy(double_data[is.fitness == TRUE & is.reads0 == TRUE])
  }
  
  ### calculate a first approximation of the fitness surface via loess
  # use at max 10^5 variants, otherwise this can be very slow
  subDT = DT[sample(nrow(DT),min(c(100000,nrow(DT)))),.(fitness1,fitness2,fitness)]
  F_fit_loess_model = loess(fitness ~ fitness1 + fitness2,data=subDT,span=0.2)
  # extrapolate loess fit to all variants
  double_data[is.fitness == TRUE & is.reads0 == TRUE,
              F_fit_loess := predict(F_fit_loess_model,newdata = .SD),,.SDcols = c("fitness1","fitness2")]
  
  # predict loess surface fitness for all variants in DT data.table 
  DT[,F_fit_loess := predict(F_fit_loess_model,newdata = .SD),,.SDcols = c("fitness1","fitness2")]
  # calculate a loess-corrected fitness
  DT[,fitness_norm := fitness-F_fit_loess]
  
  ### correct the loess approximated fitness surface by quantile surfaces
  Nq = 100 # grid points along each axis
  Nv = 500000 # max # of variants used in estimation of surface
  span = max(c(0.01,500/nrow(DT))) # fraction of nearest neighbours to use for median calculation
  
  ## calculate quantile fitness surfaces
  List = quantile_fitness_surface_adaptive(DT,Nq,Nv,span,Q)
  
  double_data[is.fitness == TRUE & is.reads0==TRUE,
          F_fit_median := predict(List$F_median_fit,newdata = .SD) + F_fit_loess,,.SDcols = c("fitness1","fitness2")]
  double_data[is.fitness == TRUE & is.reads0==TRUE,
          F_fit_lower := predict(List$F_lower_fit,newdata = .SD) + F_fit_loess,,.SDcols = c("fitness1","fitness2")]
  double_data[is.fitness == TRUE & is.reads0==TRUE,
          F_fit_upper := predict(List$F_upper_fit,newdata = .SD) + F_fit_loess,,.SDcols = c("fitness1","fitness2")]
  
  
  ### calculate error of (quantitative) epistasis estimate (double mutant fitness - median fitness surface)
  if (epistasis_error_from_slopes) {
    #calculate slope of median surface to estimate error propagation from singles
    f1 = predict(List$F_median_fit,newdata = double_data[is.fitness == TRUE & is.reads0==TRUE,.(fitness1 = fitness1 + 0.01,fitness2)]) +
      predict(F_fit_loess_model,newdata = double_data[is.fitness == TRUE & is.reads0==TRUE,.(fitness1 = fitness1 + 0.01,fitness2)])
    double_data[is.fitness == TRUE & is.reads0==TRUE,slope1 := abs(F_fit_median - f1)/0.01]
    
    f2 = predict(List$F_median_fit,newdata = double_data[is.fitness == TRUE & is.reads0==TRUE,.(fitness1,fitness2 = fitness2+ 0.01)]) +
      predict(F_fit_loess_model,newdata = double_data[is.fitness == TRUE & is.reads0==TRUE,.(fitness1,fitness2 = fitness2 + 0.01)])
    double_data[is.fitness == TRUE & is.reads0==TRUE,slope2 := abs(F_fit_median - f2)/0.01]
    
    #from this calculate epistasis error via error propagation
    double_data[,sigmaE := sqrt(sigma^2 + slope1^2 * sigma1^2 + slope2^2 * sigma2^2)]
  } else {
    #otherwise, just add variances of single and double mutant fitness estimates
    double_data[,sigmaE := sqrt(sigma^2 + sigma1^2 + sigma2^2)]
  }
  
  
  ########################################################
  ################# plot fitness surfaces ################
  ########################################################
  plot_fitness_surface(double_data,F_fit_loess_model,List,dataset_dir,prefix)
  
  
  #####################################################################
  ######## define data subsets for positive/negative epistasis ########
  #####################################################################
  
  ## estimate width (/95percent quantile) of background fitness distribution (due to measurement limit of fitness assay)
  # take only double mutants with expected fitness significantly below lower fitness bound 
  background_cutoff = double_data[is.fitness==TRUE & is.reads0 == TRUE & 
                              (fitness1 + xsig*sigma1 + fitness2 + xsig*sigma2) < lower_bound_F,
                              quantile(fitness,probs=0.95,na.rm=T)]
  
  ## mark variants for positive epistasis analysis
  double_data[is.fitness==TRUE & is.reads0 == TRUE,pos_epistasis := FALSE]
  double_data[is.fitness==TRUE & is.reads0 == TRUE &
            F_fit_upper < 0 & #q95 surface smaller than wild-type fitness
            (fitness1 + xsig*sigma1 < 0 | fitness2 + xsig*sigma2 < 0) & # f1 or f2 below wild-type fitness
            !(fitness1 - xsig*sigma1 + fitness2 - xsig*sigma2 > 0), #expected fitness not above wild-type fitness
          pos_epistasis := TRUE]
  
  ## mark variants for negative epistasis analysis
  double_data[is.fitness==TRUE & is.reads0 == TRUE,neg_epistasis := FALSE]
  double_data[is.fitness==TRUE & is.reads0 == TRUE &
            F_fit_lower > background_cutoff & # q5 surface larger than q95 of background
            fitness1 - xsig*sigma1 > lower_bound_F &
            fitness2 - xsig*sigma2 > lower_bound_F & # f1&f2 above lower fitness limit
            !(fitness1 - xsig*sigma1 + fitness2 - xsig*sigma2 > 0), #expected fitness not above wild-type fitness
          neg_epistasis := TRUE]
  
  
  ##### write double_data data.table to file
  write.table(x = double_data, file = paste0(dataset_dir,"processed_data/",output_filename),
              quote = F,row.names = F, col.names = T)
  
  #and return it
  return(double_data)
}

### calculate quantile fitness surfaces
quantile_fitness_surface_adaptive = function(DT,Nq,Nv,span,Q) {
  
  # calculate quantile surface approximation on regular grid (defined by quantiles) given a sampled number of variants
  
  #define grid vector
  q = seq(min(c(quantile(DT[,fitness1],probs=0,na.rm=T),quantile(DT[,fitness2],probs=0,na.rm=T))),
          max(c(quantile(DT[,fitness1],probs = 1,na.rm=T),quantile(DT[,fitness2],probs = 1,na.rm=T)))+0.015,
                length.out=Nq)
  
  #Fq: data.table for surface values at each gridpoint
  Fq = data.table(fitness1=rep(q,length(q)),fitness2=rep(q,each=length(q)))
  
  #initialize different surface columns (median, upper, lower)
  Fq[,F_median := as.numeric(NA)]
  Fq[,F_lower := as.numeric(NA)]
  Fq[,F_upper := as.numeric(NA)]
  
  #downsample variants to value given by Nv (or use all variants if number of variants if smaller Nv)
  subDT = DT[sample(nrow(DT),min(c(Nv,nrow(DT)))),.(fitness1,fitness2,fitness_norm)]
  
  # run this in parallel for all grid points
  require(parallel)
  # Use the detectCores() function to find the number of cores in system
  no_cores <- detectCores()-1
  clust <- makeCluster(no_cores) 
  # make variables available to each core's workspace
  clusterExport(clust, list("subDT","span","Fq","Q"),envir = environment())
  helper = parSapply(clust,X = 1:Nq^2, surface_at_gridpoint)
  stopCluster(clust)
  #transfer results from helper to Fq
  Fq[,F_median := helper[1,]]
  Fq[,F_lower := helper[2,]]
  Fq[,F_upper := helper[3,]]
  
  # > loess fit regular grid surface
  List = list()
  List$F_median_fit = loess(F_median ~ fitness1 + fitness2,data = Fq,span=0.2)
  List$F_lower_fit = loess(F_lower ~ fitness1 + fitness2,data = Fq,span=0.2)
  List$F_upper_fit = loess(F_upper ~ fitness1 + fitness2,data = Fq,span=0.2)
  
  List$Fq = Fq
  return(List)
}

#calculate fitness quantiles for nearest neighbours of gridpoint
surface_at_gridpoint = function(i) {
  A = unlist(subDT[,.(D=sqrt((fitness1-Fq[i,fitness1])^2 +(fitness2-Fq[i,fitness2])^2),fitness_norm)][D < quantile(D,probs=span,na.rm=T),
                                                                                               .(quantile(fitness_norm,p=c(0.5,Q,1-Q),na.rm=T))])
  return(A)
}

#################### plot surfaces #################### 
plot_fitness_surface = function(double_data,F_fit_loess_model,List,dataset_dir,prefix) {
  
  #range of data to plot, omitting long tails in single mutant space
  xyrange = double_data[is.fitness == T & is.reads0==T,quantile(fitness1,c(0.005,0.995),na.rm=T)]
  xy = seq(xyrange[1],
           xyrange[2],
           abs(diff(xyrange))/25)

  Fd_pred = predict(List$F_median_fit,data.frame(fitness1=rep(xy,length(xy)),fitness2=rep(xy,each=length(xy))))
  Fd_pred2 =  predict(F_fit_loess_model,data.frame(fitness1=rep(xy,length(xy)),fitness2=rep(xy,each=length(xy))))
  
  Fd_pred_q05 = predict(List$F_lower_fit,data.frame(fitness1=rep(xy,length(xy)),fitness2=rep(xy,each=length(xy))))
  Fd_pred_q95 = predict(List$F_upper_fit,data.frame(fitness1=rep(xy,length(xy)),fitness2=rep(xy,each=length(xy))))
  
  xyz = double_data[is.fitness == T & is.reads0==T,between(fitness1,xyrange[1],xyrange[2]) & between(fitness2,xyrange[1],xyrange[2]),
                    .(x=fitness1,y=fitness2,z=fitness,below_q05=fitness < F_fit_lower,above_q95 = fitness>F_fit_upper)]
  
  #number of points to plot, 10k is sufficient, otherwise PDF becomes very large
  r = sample(x = nrow(xyz),size = min(c(10000,nrow(xyz))))
  x=xyz$x
  y=xyz$y
  z=xyz$z
  
  #axis limits, adjust
  xlim_plot = ylim_plot = c(xyrange[1] - 0.1*diff(xyrange),xyrange[2] + 0.1*diff(xyrange))
  zlim_plot = quantile(xyz$z,probs = c(0.005,0.995),na.rm = T)
  #plot angles
  theta_plot = c(15,55)
  phi_plot = 15 
  for (idx in 1:length(theta_plot)) {
    #### upper and lower surface with points
    pdf(paste0(dataset_dir, "results/epistasis/",prefix,"epistasis_surface",idx,".pdf"), useDingbats=FALSE)
    a=persp(xy,xy,matrix(Fd_pred2+Fd_pred_q05,nrow=length(xy),ncol=length(xy)),
            xlab="single mutant fitness 1",ylab="single mutant fitness 2",zlab="double mutant fitness",
            xlim = xlim_plot,ylim = ylim_plot,zlim = zlim_plot, 
            theta = theta_plot[idx], phi = phi_plot, 
            col=NA, nticks=5,ticktype="detailed",expand=0.75, box=TRUE)
    b=trans3d(xyz[intersect(r,which(below_q05==T))]$x,
              xyz[intersect(r,which(below_q05==T))]$y,
              xyz[intersect(r,which(below_q05==T))]$z,a)
    points(b$x,b$y,col=rgb(1,0.1,0.1),pch=16,cex=0.75)
    par(new=TRUE)
    a=persp(xy,xy,matrix(Fd_pred2+Fd_pred_q05,nrow=length(xy),ncol=length(xy)),
            xlab="single mutant fitness 1",ylab="single mutant fitness 2",zlab="double mutant fitness",
            xlim = xlim_plot,ylim = ylim_plot,zlim = zlim_plot,
            theta = theta_plot[idx], phi = phi_plot,
            col=NA, nticks=5,ticktype="detailed",expand=0.75, box=TRUE)
    par(new=TRUE)
    b=trans3d(xyz[intersect(r,which(below_q05==F & above_q95==F))]$x,
              xyz[intersect(r,which(below_q05==F & above_q95==F))]$y,
              xyz[intersect(r,which(below_q05==F & above_q95==F))]$z,a)
    points(b$x,b$y,col=rgb(1,0.5,0.5),pch=16,cex=0.75)
    par(new=TRUE)
    a2=persp(xy,xy,matrix(Fd_pred2+Fd_pred_q95,nrow=length(xy),ncol=length(xy)),
             xlab="single mutant fitness 1",ylab="single mutant fitness 2",zlab="double mutant fitness",
             xlim = xlim_plot,ylim = ylim_plot,zlim = zlim_plot,
             theta = theta_plot[idx], phi = phi_plot,
             col=NA, nticks=5,ticktype="detailed",expand=0.75, box=TRUE)
    par(new=TRUE)
    b=trans3d(xyz[intersect(r,which(above_q95==T))]$x,
              xyz[intersect(r,which(above_q95==T))]$y,
              xyz[intersect(r,which(above_q95==T))]$z,a)
    points(b$x,b$y,col=rgb(0.5,0.9,0.1),pch=16,cex=0.75)
    dev.off()
  }
}










