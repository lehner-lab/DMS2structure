
###################################################################
##### score kernel at every position along diagonal of PWI data ###
###################################################################

kernel_structure_propensity <- function(PWI, 
  kernel, 
  dataset_dir,
  prefix = "",
  Nsamples = 10000, 
  debug_this = F, 
  rand_strategy = c("all_data", "within_kernel", "kernal_width")) {
  
  ### variables 
  # PWI: pairwise interaction score data.table; except for Pos1 and Pos2 this should only contain the scores that SS elements should be predicted from
  # kernel: pairwise interaction score matrix; symmetric numeric matrix with NAs on the diagonal
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/processed_data/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # Nsamples: number of randomized controls to compare SS propensity against
  # debug_this: if TRUE, function will stop at certain points in scripts in order to understand bugs
  
  require(data.table)
  require(metap)
  
  #Initialize result list
  ss_data_list = list()

  #Which scores should be used for prediction?
  eval_cols = setdiff(names(PWI),c("Pos1","Pos2","WT_AA1","WT_AA2","NposE","NnegE"))

  for (eval_cols_idx in seq_along(eval_cols)) {
    print(eval_cols[eval_cols_idx])
    ss_data = copy(PWI[Pos1<=Pos2,.(Pos1,Pos2,input = .SD),,.SDcols = eval_cols[eval_cols_idx]])  
    setkey(ss_data,Pos1,Pos2)
    #Perpendicular distance from diagonal
    ss_data[,pos_perp := abs(Pos1-Pos2)]
    #position range for prediction
    kernel_length <- dim(kernel)[1]
    data_range = c(min(c(ss_data$Pos1,ss_data$Pos2)),max(c(ss_data$Pos1,ss_data$Pos2)))
    data_length = length(data_range[1]:data_range[2])
    pos_range = c(min(c(ss_data$Pos1,ss_data$Pos2))-(kernel_length-2),max(c(ss_data$Pos1,ss_data$Pos2))-1)

    if (debug_this) {browser()}
    
    #Kernel structure propensities
    set.seed(1603)
    for (i in pos_range[1]:pos_range[2]) {

      #Square distances from center position
      ss_data[,within_kernel := (Pos1-i)<kernel_length & (Pos1-i)>=0 & (Pos2-i)<kernel_length & (Pos2-i)>=0]

      #Construct kernel matrix
      kernel_mat <- matrix(NA, nrow=data_length, ncol=data_length)
      colnames(kernel_mat) <- data_range[1]:data_range[2]
      rownames(kernel_mat) <- colnames(kernel_mat)
      within_kernel_mat <- kernel_mat
      within_kernel_mat[is.na(within_kernel_mat)] <- FALSE
      for(j in 1:kernel_length){
        for(k in 1:kernel_length){
          j_shift <- j + i - 1
          k_shift <- k + i - 1
          if(j_shift>=data_range[1] & j_shift<=data_range[2] & k_shift>=data_range[1] & k_shift<=data_range[2]){
            kernel_mat[j_shift, k_shift] <- kernel[j, k]
            if(j_shift!=k_shift){
              within_kernel_mat[j_shift, k_shift] <- TRUE
            }
          }
        }
      }

      #Determine kernel weights
      ss_data[Pos1 %in% data_range[1]:data_range[2] & Pos2 %in% data_range[1]:data_range[2],within_data := T]
      ss_data[!is.na(ss_data$within_data), kernel_weight := kernel_mat[cbind(Pos1,Pos2)]]
      ss_data[!is.na(ss_data$within_data), within_kernel := within_kernel_mat[cbind(Pos1,Pos2)]]
      
      #Calculate kernel smoothed value for true data
      if(dim(ss_data[Pos1==i & Pos2==i])[1]==0){
        ss_data <- rbind(ss_data, list("Pos1"=i, "Pos2"=i, "within_kernel"=F, "kernel_score" = ss_data[within_kernel==T,sum(input*kernel_weight,na.rm=T)]), fill = T)
      }else{
        ss_data[Pos1==i & Pos2==i,kernel_score := ss_data[within_kernel==T,sum(input*kernel_weight,na.rm=T)]]
      }
      
      #Calculate kernel smoothed value for random distributions
      B = copy(ss_data[within_kernel==T,.(kernel_weight,input)])
      if(rand_strategy=="all_data"){
        sample_matrix = matrix(sample(ss_data[Pos1!=Pos2,c(input)],(nrow(B))*Nsamples,replace = T),nrow = nrow(B),ncol=Nsamples)
      }
      if(rand_strategy=="kernal_width"){
        sample_matrix = matrix(sample(ss_data[Pos1!=Pos2 & pos_perp<=max(ss_data[within_kernel==T,pos_perp]),c(input)],(nrow(B))*Nsamples,replace = T),nrow = nrow(B),ncol=Nsamples)
      }
      if(rand_strategy=="within_kernel"){
        sample_matrix = matrix(sample(ss_data[Pos1!=Pos2 & within_kernel==T,c(input)],(nrow(B))*Nsamples,replace = T),nrow = nrow(B),ncol=Nsamples)
      }
      kernel_sampled = colSums(sample_matrix * matrix(rep(t(B[,kernel_weight]),Nsamples),nrow=nrow(B),ncol=Nsamples),na.rm=T)
      
      #P-value for true value
      ss_data[Pos1==i & Pos2==i,kernel_p := sum(kernel_sampled >= kernel_score)/Nsamples]      
    }
    #Avoid -Inf if logging p values by setting those positions smaller than all random samples to smallest non-zero pvalue
    ss_data[kernel_p == 0 ,kernel_p := 1/Nsamples]
    #Save
    ss_data[,(paste0(eval_cols[eval_cols_idx], "_kernel_score")) := kernel_score]
    ss_data[,(paste0(eval_cols[eval_cols_idx], "_kernel_p")) := kernel_p]
    ss_data[,Pos := Pos1]
    #Restrict to desired rows and columns, sort, save
    ss_data <- ss_data[Pos1==Pos2,.SD,.SDcols=names(ss_data)[grep("^Pos$|_kernel_score$|_kernel_p$", names(ss_data))]]
    setkey(ss_data, Pos)
    ss_data_list[[eval_cols[eval_cols_idx]]] <- ss_data
  }

  #Merge DT lists
  ss_data_merge <- Reduce(function(...) merge(..., all = T), ss_data_list)
  #Write to file
  write.table(paste0(dataset_dir,"processed_data/",prefix,"kernel_structure_propensity.txt"),
              x = ss_data_merge,quote = F,row.names = F,col.names = T)
  return(ss_data_merge)
}