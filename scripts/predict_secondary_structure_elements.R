######################################################################
##### predict alpha helices and beta strands with kernel smoothing ###
######################################################################
predict_secondary_structure_elements = function(PWI,
                                           dataset_dir,
                                           prefix = "",
                                           known_SS = c(),
                                           scale_long = 1/4^2,
                                           seed_size=3,
                                           p_detection_threshold = 0.05,
                                           Nsamples = 10000,
                                           debug_this = F,
                                           return_list = F) {
  
  ### variables 
  # PWI: pairwise interaction score data.table; except for Pos1 and Pos2 this should only contain the scores that SS elements should be predicted from
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # known_SS (optional): filepointer to a file with known secondary structure elements (from a PDB file) [to plot as comparision], a table with a position and a SS classifier column
  # scale_long: length scale for gaussian smoothing kernel
  # seed_size: number of positions the SSpropensities are initially aggregated over; must be an odd value !!!
  # p_detection_threshold: p-value threshold for calling a SS element
  # Nsamples: number of randomized controls to compare SS propensity against
  # debug_this: if TRUE, function will stop at certain points in scripts in order to understand bugs
  # return_list: if TRUE, ggplot2 and ss_data objects returned in addition to predicted secondary_structure (in named list)

  
  require(data.table)
  require(ggplot2)
  require(metap)
  
  #which scores should be used for prediction?
  eval_cols = setdiff(names(PWI),c("Pos1","Pos2","WT_AA1","WT_AA2","NposE","NnegE"))

  #Initialise list of returned data
  saved_objects <- list(
    "plot_objects" = list(), 
    "secondary_structure_score" = list())

  for (eval_cols_idx in seq_along(eval_cols)) {
    print(eval_cols[eval_cols_idx])
    ss_data = copy(PWI[Pos1<=Pos2,.(Pos1,Pos2,input = .SD),,.SDcols = eval_cols[eval_cols_idx]])  
    setkey(ss_data,Pos1,Pos2)
    #position range for prediction
    pos_range = c(min(c(ss_data$Pos1,ss_data$Pos2)),max(c(ss_data$Pos1,ss_data$Pos2)))
    
    #compute diagonal/perpendicular coordinates
    ss_data[,pos_diag := (Pos1+Pos2)/2]
    ss_data[,pos_perp := abs(Pos1-Pos2)/2] #this is half the actual distance between positions in a pair (for consistency with pos_diag)
    if (debug_this) {browser()}
    
    ### secondary structure propensities
    set.seed(1603)
    for (i in pos_range[1]:pos_range[2]) {
      j = i
      
      # hamming distances from center position
      ss_data[,ham := abs(Pos1-i) + abs(Pos2-j)]
      ss_data[,ham_perp := abs(Pos1-i - (Pos2-j))]
      ss_data[,ham_diag := abs(Pos1-i + Pos2-j)]
      
      ####################################
      ####### detect alpha helices ####### 
      ####################################
      
      # compute kernel weights
      if (i > pos_range[1]) {ss_data[,alpha_weight := NULL]}
      ss_data[ham_perp <= 5 & ham_diag < 12,alpha_weight := (cos(ham_perp*2*pi/3.6)+1/3) * exp(-scale_long*ham_diag^2)]
      ss_data[Pos1==Pos2,alpha_weight := NA]
      ss_data[is.na(input),alpha_weight := NA]
      
      # calculate kernel smoothed value for true data
      ss_data[Pos1==i & Pos2==i,alpha_score := ss_data[ham_perp <= 5 & ham_diag < 12,sum(input*alpha_weight,na.rm=T)]]
      
      # calculate kernel smoothed value for random distributions
      B = copy(ss_data[ham_perp <= 5 & ham_diag <= 12,.(ham,alpha_weight,input)])
      setkey(B,ham)
      sample_matrix = matrix(sample(ss_data[Pos1!=Pos2,c(input)],(nrow(B))*Nsamples,replace = T),nrow = nrow(B),ncol=Nsamples)
      alpha_sampled = colSums(sample_matrix * matrix(rep(t(B[,alpha_weight]),Nsamples),nrow=nrow(B),ncol=Nsamples),na.rm=T)
      
      # p value for true value
      ss_data[Pos1==i & Pos2==i,alpha_p := sum(alpha_sampled >= alpha_score)/Nsamples]
      
      ####################################
      ####### detect beta strands  ####### 
      ####################################
      
      # compute kernel weights
      if (i > pos_range[1]) {ss_data[,beta_weight := NULL]}
      ss_data[ham_perp <= 2,beta_weight := ((ham_perp+1) %% 2 - 1/3)*exp(-scale_long*ham_diag^2)]
      ss_data[ham_perp == 0,beta_weight := beta_weight * 2]
      ss_data[Pos1==Pos2,beta_weight := NA]
      ss_data[is.na(input),beta_weight := NA]
      
      # calculate kernel smoothed value for true data
      ss_data[Pos1==i & Pos2==i,beta_score := ss_data[ham_perp <= 2 & ham_diag <= 12,sum(input*beta_weight,na.rm=T)]]
      
      # calculate kernel smoothed value for random distributions
      B = copy(ss_data[ham_perp <= 2 & ham_diag <= 12,.(ham,beta_weight,input)])
      setkey(B,ham)
      sample_matrix = matrix(sample(ss_data[Pos1!=Pos2,c(input)],(nrow(B))*Nsamples,replace = T),nrow = nrow(B),ncol=Nsamples)
      beta_sampled = colSums(sample_matrix * matrix(rep(t(B[,beta_weight]),Nsamples),nrow=nrow(B),ncol=Nsamples),na.rm=T)
      
      # p value for true value
      ss_data[Pos1==i & Pos2==i,beta_p := sum(beta_sampled >= beta_score)/Nsamples]
      
    }
    
    #avoid -Inf if logging p values by setting those positions smaller than all random samples to smallest non-zero pvalue
    ss_data[alpha_p == 0 ,alpha_p := 1/Nsamples]
    ss_data[beta_p == 0 ,beta_p := 1/Nsamples]
    
    ###########################################################
    ### call secondary structure elements from propensities ###
    ###########################################################
    
    ### get alpha helices and beta strand p values
    setkey(ss_data,Pos1)
    ss_strands = ss_data[Pos1==Pos2,.(Pos1,alpha_p,beta_p)]
    
    # compute sumlog (combined p-values) for seeds of alpha helix
    if (seed_size > 1) {
      for (mid_idx in ss_strands[!is.na(alpha_p),Pos1]) {
        ss_strands[Pos1==mid_idx,alpha_p_seed := ss_strands[between(Pos1,mid_idx-(seed_size-1)/2,mid_idx+(seed_size-1)/2),
                                                            ifelse(sum(!is.na(alpha_p)) > 1,sumlog(alpha_p[!is.na(alpha_p)])$p,alpha_p[!is.na(alpha_p)])]]  
      }
    } else {
      ss_strands[,alpha_p_seed := alpha_p]
    }
    
    # compute sumlog for seeds of beta strands
    if (seed_size > 1) {
      for (mid_idx in ss_strands[!is.na(beta_p),Pos1]) {
        ss_strands[Pos1==mid_idx,beta_p_seed := ss_strands[between(Pos1,mid_idx-(seed_size-1)/2,mid_idx+(seed_size-1)/2),
                                                           ifelse(sum(!is.na(beta_p)) > 1,sumlog(beta_p[!is.na(beta_p)])$p,beta_p[!is.na(beta_p)])]]  
      }
    } else {
      ss_strands[,beta_p_seed := beta_p]
    }
    
    #set p-values NA if other structure is more probable
    ss_strands[,beta_p0 := beta_p]
    ss_strands[,alpha_p0 := alpha_p]
    ss_strands[,beta_p_seed0 := beta_p_seed]
    ss_strands[,alpha_p_seed0 := alpha_p_seed]
    if (seed_size > 1) {
      if (ss_strands[1+(seed_size-1)/2]$alpha_p_seed < ss_strands[1+(seed_size-1)/2]$beta_p_seed) {
        ss_strands[1,':=' (beta_p = NA, beta_p_seed = NA)] } else {
          ss_strands[1,':=' (alpha_p = NA, alpha_p_seed = NA)] }
      if (ss_strands[.N-(seed_size-1)/2]$alpha_p_seed < ss_strands[.N-(seed_size-1)/2]$beta_p_seed) {
        ss_strands[.N,':=' (beta_p = NA, beta_p_seed = NA)] } else {
          ss_strands[.N,':=' (alpha_p = NA, alpha_p_seed = NA)] }
    }
    ss_strands[alpha_p_seed < beta_p_seed,':=' (beta_p_seed = NA, beta_p = NA)]
    ss_strands[beta_p_seed < alpha_p_seed,':=' (alpha_p_seed = NA, alpha_p = NA)]
    
    #delete stretches smaller 5 for beta strands, set NA
    ss_strands[ss_strands[,.(stretch = rleid(!is.na(alpha_p_seed)),not_na = !is.na(alpha_p_seed))][,.(short = .N<5 & not_na == T),stretch][,short],
               ':=' (alpha_p =NA,alpha_p_seed = NA)]
    #delete stretches smaller 3 for beta strands, set NA
    ss_strands[ss_strands[,.(stretch = rleid(!is.na(beta_p_seed)),not_na = !is.na(beta_p_seed))][,.(short = .N<3 & not_na == T),stretch][,short],
               ':=' (beta_p = NA,beta_p_seed = NA)]
    
    
    
    ## alpha helices:
    # identify most significant stretches from seeds
    helper = identify_expand_seeds(ss_strands[,.(pos=Pos1,p_ind=alpha_p,p_seed=alpha_p_seed)],seed_size,p_detection_threshold)
    # merge
    ss_strands = merge(ss_strands,helper[,.(Pos1=pos,alpha_strand=strand,alpha_strand_p=p_strand)],by="Pos1")
    
    ## beta strands:
    # identify most significant stretches from seeds
    helper = identify_expand_seeds(ss_strands[,.(pos=Pos1,p_ind=beta_p,p_seed = beta_p_seed)],seed_size,p_detection_threshold)
    #merge
    ss_strands = merge(ss_strands,helper[,.(Pos1=pos,beta_strand=strand,beta_strand_p=p_strand)],by="Pos1")
    #Save
    if(return_list){
      saved_objects[["secondary_structure_score"]][[eval_cols[eval_cols_idx]]] <- copy(ss_strands)
    }
    
    ### record predictions across input data
    if (eval_cols_idx == 1) {
      if (debug_this) {browser()}
      secondary_structure = ss_strands[,.(Pos = Pos1,ss = ifelse(is.na(alpha_strand) & is.na(beta_strand),"C",ifelse(!is.na(alpha_strand),"H","E")))]
      names(secondary_structure)[2] = eval_cols[1]
    } else {
      secondary_structure = merge(secondary_structure,
                                  ss_strands[,.(Pos = Pos1,ss = ifelse(is.na(alpha_strand) & is.na(beta_strand),"C",ifelse(!is.na(alpha_strand),"H","E")))],
                                  by = "Pos")
      names(secondary_structure)[eval_cols_idx+1] = eval_cols[eval_cols_idx]
    }
    
    
    #############
    ## compare to known_SS
    #############
    if (!is.null(known_SS)) {
      known_ss_DT = fread(known_SS)
      names(known_ss_DT)[2] = "SS"
      known_ss_DT[,rleidx := rleid(SS)]
      known_ss_DT = rbind(known_ss_DT,data.table(Pos = nrow(known_ss_DT)+1, SS = "C", rleidx = max(known_ss_DT$rleidx)+1))
      known_ss_DT = rbind(known_ss_DT,data.table(Pos = nrow(known_ss_DT)+1, SS = "E", rleidx = max(known_ss_DT$rleidx)+1))
      known_ss_DT = rbind(known_ss_DT,data.table(Pos = nrow(known_ss_DT)+1, SS = "H", rleidx = max(known_ss_DT$rleidx)+1))
    }
    
    if (debug_this) {browser()}
    ##################################################
    ##### plot secondary structure element predictions
    require(cowplot)
    theme_set(theme_classic())
    P1 = ggplot(data=ss_strands) + 
      geom_line(aes(Pos1,y=alpha_p_seed0),color="darkgreen",linetype=2) +
      geom_line(aes(Pos1,y=alpha_p_seed),color="darkgreen") +
      
      geom_line(aes(Pos1,y=beta_p_seed0),color="red",linetype=2) +
      geom_line(aes(Pos1,y=beta_p_seed),color="red") +
      
      geom_hline(yintercept=0.05,linetype=3) +
      geom_rect(data=unique(ss_strands[!is.na(alpha_strand),.(xmin=min(Pos1)-0.5,xmax=max(Pos1)+0.5,ymin=max(alpha_strand_p,10^-4),ymax=1),alpha_strand]),
                inherit.aes = F,
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,group=alpha_strand),fill="green",alpha=0.2) +
      geom_rect(data=unique(ss_strands[!is.na(beta_strand),.(xmin=min(Pos1)-0.5,xmax=max(Pos1)+0.5,ymin=max(beta_strand_p,10^-4),ymax=1),beta_strand]),
                inherit.aes = F,
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,group=beta_strand),fill="orange",alpha=0.2) +
      scale_x_continuous(breaks=seq(5,pos_range[2],5),expand = c(0,0)) +
      coord_cartesian(xlim = c(pos_range[1]-0.5,pos_range[2]+0.5)) +
      # scale_y_log10(breaks = c(10^-seq(-10,0))) +
      scale_y_log10(breaks = c(10^-seq(0,ss_strands[,-log10(min(c(alpha_p_seed,beta_p_seed,10^-4),na.rm=T))])),expand = c(0.01,0)) +
      labs(y="p value",title = eval_cols[eval_cols_idx])
    
    if (!is.null(known_SS)) { #add secondary structure
      P1 = P1 +  
        geom_segment(data = known_ss_DT[,.(start = min(Pos)-0.5,end = max(Pos)+0.5,ss=unique(SS)),rleidx],
                     aes(x=start,y=1.5,xend=end,yend=1.5,color=ss,size=ss),show.legend = F) +
        scale_size_manual(breaks = c("C","E","H"),values = c(0.5,1.5,1.5)) +
        scale_color_manual(breaks = c("C","E","H"),values = c("black","red","darkgreen"))
    }
    
    ##################################################
    ###### plot smoothed data around diagonal
    setkey(ss_data,Pos1,Pos2)
    
    for (i in pos_range[1]:pos_range[2]) {
      for (j in i:pos_range[2]) {
        ss_data[,ham := abs(Pos1-i) + abs(Pos2-j)]
        ss_data[,ham_perp := abs(Pos1-i - (Pos2-j))]
        ss_data[,ham_diag := abs(Pos1-i + Pos2-j)]
        
        ss_data[,weight := as.double(NA)]
        ss_data[ham_perp==0,weight := exp(-scale_long*ham_diag^2)]
        ss_data[Pos1==Pos2,weight := NA]
        ss_data[is.na(input),weight := NA]
        #calculate true value
        ss_data[Pos1==i & Pos2==j,score := ss_data[ham_perp==0,sum(input*weight,na.rm=T)/sum(weight,na.rm=T)]]
      }
    }
    #average
    avg = 1.5
    for (i in ss_data[,unique(ham_diag)]) {
      ss_data[ham_diag == i & ham_perp < 8,score_norm := score-ss_data[between(ham_diag,i-avg,i+avg) & ham_perp < 10,mean(score,na.rm=T)]]
      ss_data[ham_diag == i & ham_perp < 8,input_norm := input-ss_data[between(ham_diag,i-avg,i+avg) & ham_perp < 10 & abs(input) != Inf,mean(input,na.rm=T)]]
    }
    #limit data range for better comparability
    cutoff = quantile(c(ss_data[Pos1<Pos2 & !is.na(input_norm),abs(input_norm)],ss_data[!is.na(score_norm),abs(score_norm)]),0.95)
    ss_data[input_norm > cutoff, input_norm := cutoff]
    ss_data[input_norm < -cutoff, input_norm := -cutoff]
    ss_data[score_norm > cutoff, score_norm := cutoff]
    ss_data[score_norm < -cutoff, score_norm := -cutoff]
    
    P2 = ggplot(ss_data[ham_perp < 8 & ham_perp != 0]) +
      geom_raster(aes(Pos1,Pos2,fill=input_norm)) +
      geom_raster(aes(Pos2,Pos1,fill=score_norm)) +
      scale_fill_gradient2(midpoint=0,low="tomato3",high="steelblue3",na.value = "white") +
      scale_y_continuous(breaks=seq(1,8,1),expand = c(0,0)) + 
      scale_x_continuous(breaks=seq(5,pos_range[2],5),expand = c(0,0)) +
      coord_cartesian(xlim = c(pos_range[1]-0.5,pos_range[2]+0.5),ylim = c(pos_range[1]-0.5,pos_range[2]+0.5)) +
      labs(x = "position",y = "diagonal position",fill = "score") +
      geom_segment(data=unique(ss_strands[!is.na(alpha_strand),.(x=min(Pos1)-0.5,xend=max(Pos1)+0.5,y=min(Pos1)-0.5,yend=max(Pos1)+0.5),alpha_strand]),
                   inherit.aes = F,
                   aes(x=x,xend=xend,y=y,yend=yend,group=alpha_strand),color="darkgreen",size=1.5) +
      geom_segment(data=unique(ss_strands[!is.na(beta_strand),.(x=min(Pos1)-0.5,xend=max(Pos1)+0.5,y=min(Pos1)-0.5,yend=max(Pos1)+0.5),beta_strand]),
                   inherit.aes = F,
                   aes(x=x,xend=xend,y=y,yend=yend,group=beta_strand),color="red",size=1.5)
    
    plot_grid(P1,P2,nrow=1)
    
    ggsave(file=paste0(dataset_dir,"results/secondary_structure/",prefix,eval_cols[eval_cols_idx],"_SSelements.pdf"),width=8.5,height=4)

    if(return_list){
      saved_objects[["plot_objects"]][[eval_cols[eval_cols_idx]]] <- list("P1" = P1, "P2" = P2)
    }
  }
  
  write.table(paste0(dataset_dir,"processed_data/",prefix,"secondary_structure_prediction.txt"),
              x = secondary_structure,quote = F,row.names = F,col.names = T)
  
  #Return objects
  if(return_list){
    saved_objects[["secondary_structure"]] <- secondary_structure
    return(saved_objects)
  }else{
    return(secondary_structure)
  }
}

