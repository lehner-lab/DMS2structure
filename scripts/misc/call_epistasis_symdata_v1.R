
#################################################################
### a modified version of epistasis classification for "symmetrical" data with detrimental and beneficial variants ###
#################################################################

call_epistasis_symdata_v1 = function(double_data,lower_bound_F,upper_bound_F,output_dir,prefix = "",xsig=2,sym = T,Q = 0.05, epistasis_error_from_slopes = T) {
  
  
  #### version for "symmetrical" data with detrimental and beneficial variants
  # this has a modified version of epistasis classification and needs to know the lower and upper (!) bound of the fitness data
  
  #double_data: the doubles data.table
  #lower_bound_F: an estimate of the lower bound of the fitness assay
  #upper_bound_F: an estimate of the upper bound of the fitness assay
  #output_dir: base-output directory, like "GB1/", it will put results output_dir/results/epistasis/
  #prefix:to be added to results files (in case of running diff. versions of data from same dataset etc)
  #xsig: significance threshold for calling epistasis, 2 is fine
  #sym: logical for whether the position-position fitness map is symmetrical (only FALSE for protein-protein interactions!!)
  #Q: the percentile used for upper (1-Q) and lower (Q) fitness surface calculation
  #epistasis_error_from_slopes: logical, should sigmeE (error for epistasis values) be calculate by taking slopes of median fitness surface into account
  
  
  set.seed(1603)
  
  require(data.table)
  require(mgcv)
  require(caTools)
  
  if (sym == T) {
    DT = switch_double_DT(double_data[is.fitness == TRUE & is.reads0 == TRUE],
                          cols_switchdouble = list(c("fitness1","fitness2"),c("sigma1","sigma2")),
                          cols_double = c("fitness","sigma"))
  } else {
    DT = copy(double_data[is.fitness == TRUE & is.reads0 == TRUE])
  }
  
  ### use loess instead, 2d-gam is very inaccurate
  # browser()
  subDT = DT[sample(nrow(DT),min(c(100000,nrow(DT)))),.(fitness1,fitness2,fitness)]
  F_fit_loess_model = loess(fitness ~ fitness1 + fitness2,data=subDT,span=0.2)
  double_data[is.fitness == TRUE & is.reads0 == TRUE,
          F_fit_loess := predict(F_fit_loess_model,newdata = .SD),,.SDcols = c("fitness1","fitness2")]
  
  ### >> calculate A-B-AB surface as median surface of gam-corrected surface
  Nq = 100 # grid points along one axis
  Nv = 500000 # variants used in estimation of surface
  span = max(c(0.01,500/nrow(DT))) # fraction of nearest neighbours to use for median calculation
  # Q = 0.05: values for quantile calcultion: 0.05 and 0.95, 
  
  # predict loess surface fitness and correct fitness for it
  DT[,F_fit_loess := predict(F_fit_loess_model,newdata = .SD),,.SDcols = c("fitness1","fitness2")]
  DT = DT[,fitness_norm := fitness-F_fit_loess]
  
  # calculate quantile fitness surfaces
  List = quantile_fitness_surface_adaptive(DT,Nq,Nv,span,Q)

  double_data[is.fitness == TRUE & is.reads0==TRUE,
          F_fit_median := predict(List$F_median_fit,newdata = .SD) + F_fit_loess,,.SDcols = c("fitness1","fitness2")]
  double_data[is.fitness == TRUE & is.reads0==TRUE,
          F_fit_lower := predict(List$F_lower_fit,newdata = .SD) + F_fit_loess,,.SDcols = c("fitness1","fitness2")]
  double_data[is.fitness == TRUE & is.reads0==TRUE,
          F_fit_upper := predict(List$F_upper_fit,newdata = .SD) + F_fit_loess,,.SDcols = c("fitness1","fitness2")]
  
  if (epistasis_error_from_slopes) {
    #calculate slope of median surface to estimate error propagation from singles
    f1 = predict(List$F_median_fit,newdata = double_data[is.fitness == TRUE & is.reads0==TRUE,.(fitness1 = fitness1 + 0.01,fitness2)]) +
      predict(F_fit_loess_model,newdata = double_data[is.fitness == TRUE & is.reads0==TRUE,.(fitness1 = fitness1 + 0.01,fitness2)])
    double_data[is.fitness == TRUE & is.reads0==TRUE,slope1 := abs(F_fit_median - f1)/0.01]
    
    f2 = predict(List$F_median_fit,newdata = double_data[is.fitness == TRUE & is.reads0==TRUE,.(fitness1,fitness2 = fitness2+ 0.01)]) +
      predict(F_fit_loess_model,newdata = double_data[is.fitness == TRUE & is.reads0==TRUE,.(fitness1,fitness2 = fitness2 + 0.01)])
    double_data[is.fitness == TRUE & is.reads0==TRUE,slope2 := abs(F_fit_median - f2)/0.01]
    
    #from this calculate epistasis error
    double_data[,sigmaE := sqrt(sigma^2 + slope1^2 * sigma1^2 + slope2^2 * sigma2^2)]
  } else {
    double_data[,sigmaE := sqrt(sigma^2 + sigma1^2 + sigma2^2)]
  }
  
  #####################################################################
  ######## define data subsets for positive/negative epistasis ########
  #####################################################################
  
  #estimate the width (95percent quantile) of the lower limit "background"
  lowerlimit_background_cutoff=double_data[is.fitness==TRUE & is.reads0 == TRUE & 
                              (fitness1 + xsig*sigma1 + fitness2 + xsig*sigma2) < lower_bound_F,quantile(fitness,probs=0.95,na.rm=T)]
  
  #same for the upper limit background
  upperlimit_background_cutoff=double_data[is.fitness==TRUE & is.reads0 == TRUE & 
                                             (fitness1 - xsig*sigma1 + fitness2 - xsig*sigma2) > upper_bound_F,quantile(fitness,probs=0.95,na.rm=T)]
  
  # mark variants for positive epistasis analysis
  double_data[is.fitness==TRUE & is.reads0 == TRUE,pos_epistasis := FALSE]
  
  ## upper fitness restrictions for positive epistasis
  # only limitation here is that it is not too high into the upper background of fitness
  double_data[is.fitness==TRUE & is.reads0 == TRUE &
                F_fit_upper < upperlimit_background_cutoff,
          pos_epistasis := TRUE]
  
  # mark variants for negative epistasis analysis
  double_data[is.fitness==TRUE & is.reads0 == TRUE,neg_epistasis := FALSE]
  # #upper fitness restrictions for negative epistasis
  # only limitation here is that it is not too low into the lower background of fitness
  double_data[is.fitness==TRUE & is.reads0 == TRUE &
            F_fit_lower > lowerlimit_background_cutoff,
          neg_epistasis := TRUE]
  
  ########################################################
  ################# plot fitness surfaces ################
  ########################################################
  plot_fitness_surface(double_data,F_fit_loess_model,List,output_dir,prefix)
  
  #Epistasis score and significant positive/negative classifications
  double_data[is.fitness==TRUE & is.reads0 == TRUE,epistasis := fitness - F_fit_median]
  double_data[pos_epistasis == TRUE,pos_epistasis_sig := fitness-F_fit_upper > 0]
  double_data[neg_epistasis == TRUE,neg_epistasis_sig := fitness-F_fit_lower < 0]

  return(double_data)
}
