
#############################################################
##### Tau pattern-specificity score (Yanai, et al 2004)   ###
#############################################################

tau_specificity_score <- function(x, 
  min_length=3){

  ### variables 
  # x: vector of scores
  # min_length: minimum number of scores (otherwise returns NA)

  x <- x[!is.na(x)]
  if(length(x)>=min_length){
    if(max(x)==0){
      return(NA)
    }else{
      return(sum(1-x/max(x))/(length(x)-1))
    }
  }else{
    return(NA)
  }
}