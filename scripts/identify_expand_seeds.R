################################################################################################
############### function for expanding seeds and finding most significant stretches ############
################################################################################################
## subfunctions to predict_secondary_structure / sheets functions
identify_expand_seeds = function(any_strand,seed_size = 3,p_threshold = 0.05,bridge_dist = 1,max_extension = 1) {
  
  #starting with lowest p-value triplet, expand
  idx=1
  any_strand_temp = copy(any_strand)
  any_strand[,':=' (strand = as.integer(NA), p_strand = as.numeric(NA))]
  while (min(any_strand_temp$p_seed,na.rm=T) < p_threshold) {
    #start with best available triplet seed
    psum = any_strand_temp[,min(p_seed,na.rm=T)]
    positions = any_strand_temp[which.min(p_seed),seq(pos-(seed_size-1)/2,pos+(seed_size-1)/2,1)]
    
    #expand triplet seed
    positions_down = c(positions,min(positions)-1)
    positions_up = c(positions,max(positions)+1)
    extension = 0
    while (extension <= max_extension) {
      
      psum_down = any_strand_temp[pos %in% positions_down & !is.na(p_ind),if(.N>1){sumlog(p_ind)$p}else{p_ind}]
      psum_up = any_strand_temp[pos %in% positions_up & !is.na(p_ind),if(.N>1){sumlog(p_ind)$p}else{p_ind}]
      if (psum_down < min(psum,psum_up)) {
        positions = positions_down
        psum = psum_down
        extension = 0
      } else if (psum_up < min(psum,psum_down)) {
        positions = positions_up
        psum = psum_up
        extension = 0
      } else { #check whether extending further gives a more significant strand
        extension = extension+1
      }
      positions_down = c(positions_down,min(positions)-1-extension)
      positions_up = c(positions_up,max(positions)+1+extension)
    }
    
    #record beta strands
    if (idx == 1) {
      any_strand[pos %in% positions,':=' (strand = idx, p_strand = psum)]
      idx=idx+1
    } else {
      #merge with adjacent beta strands if dist smaller than bridge_dist
      closeby_strands = unique(any_strand[!is.na(strand)][,.(min_dist=min(abs(pos-positions)),strand),pos][min_dist<= bridge_dist,strand])
      if (length(closeby_strands) > 0) {
        if (length(closeby_strands) > 1) {
          any_strand[strand %in% closeby_strands,strand := min(closeby_strands)]
        }
        strand_idx = min(closeby_strands)
        any_strand[pos %in% positions,':=' (strand = strand_idx)]
        any_strand[strand == strand_idx, p_strand := sumlog(p_ind[!is.na(p_ind)])$p, strand]
      } else { # or record new beta strand
        strand_idx = idx
        any_strand[pos %in% positions,':=' (strand = strand_idx)]
        any_strand[strand == strand_idx, p_strand := sumlog(p_ind[!is.na(p_ind)])$p, strand]
        idx=idx+1
      }
    }
    
    #set p-values of 'used' positions to NA
    any_strand_temp[pos %in% positions,p_ind := NA]
    # any_strand_temp[between(pos,min(positions)-(seed_size-1)/2,max(positions)+(seed_size-1)/2),p_seed:=NA]
    any_strand_temp[pos %in% positions,p_seed:=NA]
  }
  
  #rearrange strand nrs
  strand_vec = unique(any_strand[!is.na(strand),strand])
  if (length(strand_vec) > 1) {
    for (i in seq_along(strand_vec)) {any_strand[strand==strand_vec[i],strand_new:=i]}
    any_strand[,strand:=strand_new]
    any_strand[,strand_new:=NULL]
  }
  return(any_strand)
}