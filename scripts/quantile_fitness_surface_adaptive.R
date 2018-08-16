### calculate quantile fitness surfaces
#subfunction for the call_epistasis class of functions
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