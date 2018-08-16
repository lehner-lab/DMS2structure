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












