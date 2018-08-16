############################################################
##### predict beta she et pairing with kernel smoothing #####
############################################################
predict_beta_sheets = function(PWI,
                               input_ss0,
                               dataset_dir,
                               prefix = "",
                               known_ss_file = c(),
                               known_bsi_file = c(),
                               scale_long = 1/4^2,
                               seed_size=3,
                               p_detection_threshold = 0.05,
                               debug_this = F,
                               Nsamples = 10000,
                               restricted_pairing = F) {
  
  
  
  ### variables 
  # PWI: pairwise interaction score data.table; except for Pos1 and Pos2 this should only contain the scores that SS elements should be predicted from
  # input_ss0: secondary structure input which the predictions should be based on [a table with a position column and a SS identifier column]
  #           - either from DMS data from the predict_secondary_structure_elements function
  #           - or from a PSIPRED or PDB file
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # known_SS_file (optional): filepointer to a file with known secondary structure elements (from PDB or PSIPRED) [to plot as comparision], a table with a position and a SS classifier column
  # known_bsi_file (optional): filepointer to a file with known beta sheet pairing info (from a PDB file) [to plot as comparision]
  # scale_long: length scale for gaussian smoothing kernel
  # seed_size: number of positions the sheet propensities are initially aggregated over; must be an odd value !!!
  # p_detection_threshold: p-value threshold for calling a sheet
  # Nsamples: number of randomized controls to compare sheet propensity against
  # debug_this: if TRUE, function will stop at certain points in scripts in order to understand bugs
  # restricted_pairing: if TRUE, algorithm will only evaluate potential beta sheet pairings in regions of predicted beta strands
  #                     if FALSE, algorithm will create new beta strands in regions not occupied by alpha helices if it sees fit
  
  
  require(data.table)
  require(ggplot2)
  require(metap)
  
  eval_cols = setdiff(names(PWI),c("Pos1","Pos2","WT_AA1","WT_AA2","NposE","NnegE"))
  
  for (eval_cols_idx in seq_along(eval_cols)) {
    print(eval_cols[eval_cols_idx])
    
    ss_data = copy(PWI[Pos1<=Pos2,.(Pos1,Pos2,input = .SD),,.SDcols = eval_cols[eval_cols_idx]])  
    setkey(ss_data,Pos1,Pos2)
    
    if (ncol(input_ss0) > 2){
      input_ss = input_ss0[,.(Pos,SS = .SD),,.SDcols = eval_cols[eval_cols_idx]]
    } else if (ncol(input_ss0) == 2) {
      input_ss = input_ss0[,.(Pos,SS = .SD),,.SDcols = 2]
    } else {
      print("number secondary structure inputs doesn't match the number of input features")
      input_ss = error
    }
    input_ss[,rleidx := rleid(SS)]
    
    #position range
    pos_range = c(min(c(ss_data$Pos1,ss_data$Pos2)),max(c(ss_data$Pos1,ss_data$Pos2)))
    
    #compute diagonal/perpendicular coordinates
    ss_data[,pos_diag := (Pos1+Pos2)/2]
    ss_data[,pos_perp := abs(Pos1-Pos2)/2] #this is half the actual distance between positions in a pair (for consistency with pos_diag)
    
    ss_data[Pos1==Pos2,input := NA]
    
    ### calculate sheet propensities on all off-diagonal elements of the interaction score matrix
    set.seed(1603)
    
    for (i in pos_range[1]:pos_range[2]) {
      for (j in i:pos_range[2]) {
        
        # hamming distances from center position
        ss_data[,ham := abs(Pos1-i) + abs(Pos2-j)]
        ss_data[,ham_perp := abs(Pos1-i - (Pos2-j))]
        ss_data[,ham_diag := abs(Pos1-i + Pos2-j)]
        
        ##########################################
        ####### detect parallel betasheets #######
        ##########################################
        
        # compute kernel weights
        ss_data[ham_perp <= 2,beta_par_weight := ((ham_perp+1) %% 2 - 1/3)*exp(-scale_long*ham_diag^2)]
        ss_data[ham_perp == 0,beta_par_weight := beta_par_weight * 2]
        ss_data[Pos1==Pos2,beta_par_weight := NA]
        ss_data[is.na(input),beta_par_weight := NA]
        
        # calculate kernel smoothed value for true data
        ss_data[Pos1==i & Pos2==j,beta_par_score := ss_data[ham_perp <= 2 & ham_diag <= 12,sum(input*beta_par_weight,na.rm=T)]]
        
        # calculate kernel smoothed value for random distributions
        B = copy(ss_data[ham_perp <= 2 & ham_diag <= 12,.(ham,beta_par_weight,input)])
        setkey(B,ham)
        sample_matrix = matrix(sample(ss_data[Pos1!=Pos2,c(input)],(nrow(B))*Nsamples,replace = T),nrow = nrow(B),ncol=Nsamples)
        beta_par_sampled = colSums(sample_matrix * matrix(rep(t(B[,beta_par_weight]),Nsamples),nrow=nrow(B),ncol=Nsamples),na.rm=T)
        
        # p value for true value
        ss_data[Pos1==i & Pos2==j,beta_par_p := sum(beta_par_sampled >= beta_par_score)/Nsamples]
        
        
        ###############################################
        ####### detect anti-parallel betasheets #######
        ###############################################
        
        if (j > i+1) {
          # compute kernel weights
          # ss_data[,beta_antipar_weight:=NULL]
          ss_data[ham_diag <= 2,beta_antipar_weight := ((ham_diag+1) %% 2 - 1/3)*exp(-scale_long*ham_perp^2)]
          ss_data[ham_diag == 0,beta_antipar_weight := beta_antipar_weight * 2]
          ss_data[Pos1==Pos2,beta_par_weight := NA]
          ss_data[is.na(input),beta_antipar_weight := NA]
          
          # calculate kernel smoothed value for true data
          ss_data[Pos1==i & Pos2==j,beta_antipar_score := ss_data[ham_diag <= 2 & ham_perp <= 12,sum(input*beta_antipar_weight,na.rm=T)]]
          
          # calculate kernel smoothed value for random distributions
          B = copy(ss_data[ham_diag <= 2 & ham_perp <= 12,.(ham,beta_antipar_weight,input)])
          setkey(B,ham)
          sample_matrix = matrix(sample(ss_data[Pos1!=Pos2,c(input)],(nrow(B))*Nsamples,replace = T),nrow = nrow(B),ncol=Nsamples)
          beta_antipar_sampled = colSums(sample_matrix * matrix(rep(t(B[,beta_antipar_weight]),Nsamples),nrow=nrow(B),ncol=Nsamples),na.rm=T)
          
          # p value for true value
          ss_data[Pos1==i & Pos2==j,beta_antipar_p := sum(beta_antipar_sampled >= beta_antipar_score)/Nsamples]
        }
      }
    }
    if (debug_this) {browser()}
    
    
    #avoid -Inf if logging p values by setting those positions smaller than all random samples to smallest non-zero pvalue
    ss_data[beta_par_p == 0 ,beta_par_p := 1/Nsamples]
    ss_data[beta_antipar_p == 0 ,beta_antipar_p := 1/Nsamples]
    
    ##########################################
    ####### detect beta sheet pairings #######
    ##########################################
    
    ##################################################
    calculate_seed_pval = function(ss_data,
                                   restricted_pairing,
                                   input_ss,seed_size = 3,
                                   blocked_space = c(),
                                   beta_par_distfromdiag = 4,
                                   beta_antipar_distfromdiag = 1) {
      
      #initialize variables
      ss_data[,sheet_par_nr := as.numeric(NA)]
      ss_data[,sheet_par_p := as.numeric(NA)]
      ss_data[,sheet_antipar_nr := as.numeric(NA)]
      ss_data[,sheet_antipar_p := as.numeric(NA)]
      
      directionality = c("par","antipar")
      for (ap in seq_along(directionality)) {
        if (directionality[ap] == "par") {
          ###################################
          ####### parallel beta sheets ######
          ###################################
          guide = ss_data[,.(Pos1,Pos2,pos_perp)][
            pos_perp > beta_par_distfromdiag,.(.N,pos_guide=pos_perp),pos_perp]
          
          pairing_helper = ss_data[Pos1 < Pos2 & pos_perp %in% guide$pos_guide,
                                   .(Pos1,Pos2,pos_aligned = pos_diag,pos_guide = pos_perp,p_ind = beta_par_p)]
        } else if (directionality[ap] == "antipar") {
          ###################################
          ##### anti-parallel beta sheets ###
          ###################################
          guide = ss_data[,.(Pos1,Pos2,pos_diag)][
            pos_diag > beta_antipar_distfromdiag,.(.N,pos_guide=pos_diag),pos_diag]
          
          pairing_helper = ss_data[Pos1 < Pos2 & pos_diag %in% guide$pos_guide & pos_perp > beta_antipar_distfromdiag,
                                   .(Pos1,Pos2,pos_aligned = pos_perp,pos_guide = pos_diag,p_ind = beta_antipar_p)]
        }
        
        # if (debug_this) {browser()}
        
        if (restricted_pairing) { #if only restricted pairing of strands is allowed, set all p-values of positions that do not belong to two beta strands to NA
          pairing_helper[!(Pos1 %in% input_ss[SS=="E",Pos] & Pos2 %in% input_ss[SS=="E",Pos]),p_ind := NA]
        } 
        
        if (length(blocked_space) > 0) { #if there are blocked spaces define, set all p-values of positions in those to NA
          for (i in 1:nrow(blocked_space)) {
            pairing_helper[between(Pos1,blocked_space[i,pos_min],blocked_space[i,pos_max]) | 
                             between(Pos2,blocked_space[i,pos_min],blocked_space[i,pos_max]),
                           p_ind := NA]
          }
        }
        
        #for each row of position pairs with equal pos_perp
        for (k in seq_along(guide$pos_guide)) {
          #extract line
          specific_line = pairing_helper[pos_guide == guide$pos_guide[k],.(pos_aligned,p_ind,p_seed=1)]
          
          #calculate seeds
          for (pairing_idx in specific_line[!is.na(pos_aligned),unique(pos_aligned)]) {
            p_vec = specific_line[between(pos_aligned,
                                          pairing_idx-(seed_size-1)/2,
                                          pairing_idx+(seed_size-1)/2),p_ind]
            if (length(p_vec)==seed_size) {
              if (seed_size > 1) {
                # specific_line[pos_aligned == pairing_idx,p_seed := sumlog(p_vec)$p]
                # specific_line[pos_aligned == pairing_idx,p_seed := ifelse(sum(!is.na(p_vec)) > 1,
                # sumlog(p_vec[!is.na(p_vec)])$p,
                # NA)]
                specific_line[pos_aligned == pairing_idx,p_seed := ifelse(sum(!is.na(p_vec)) == seed_size,
                                                                          sumlog(p_vec)$p,
                                                                          NA)]
              } else {
                specific_line[pos_aligned == pairing_idx,p_seed := p_vec]
              }
            }
          }
          
          #extend seeds
          helper = identify_expand_seeds(specific_line[,.(pos=pos_aligned,p_ind,p_seed)],seed_size)
          
          #collect info
          if (directionality[ap] == "par") {
            ss_data[Pos1<Pos2 &
                      pos_perp == guide$pos_guide[k] &
                      pos_diag %in% helper[!is.na(p_strand),pos],
                    sheet_par_p:=helper[pos==pos_diag,p_strand],pos_diag]
            max_strand = max(c(ss_data$sheet_par_nr,0),na.rm=T)
            ss_data[Pos1<Pos2 &
                      pos_perp == guide$pos_guide[k] &
                      pos_diag %in% helper[!is.na(p_strand),pos],
                    sheet_par_nr:=helper[pos==pos_diag,strand + max_strand],pos_diag]
          } else if (directionality[ap] == "antipar") {
            # check whether there's adjacent/overlapping strands reported from other sheetpairing
            for (s in helper[!is.na(strand),unique(strand)]) {
              if (nrow(ss_data[Pos1<Pos2 & pos_diag == guide$pos_guide[k] & !is.na(sheet_antipar_nr)]) > 0) {
                closeby_strands = unique(ss_data[Pos1<Pos2 & pos_diag == guide$pos_guide[k] & !is.na(sheet_antipar_nr)][,
                                                                                                                        .(min_dist=min(abs(pos_perp-helper[strand==s,pos])),sheet_antipar_nr),pos_perp][min_dist <= 1,sheet_antipar_nr])
              } else { closeby_strands = c() }
              if (length(closeby_strands) > 0) { #if yes, merge
                if (length(closeby_strands) > 1) { #if there's strands on both sides, first merge those
                  ss_data[sheet_antipar_nr %in% closeby_strands,sheet_antipar_nr := min(closeby_strands)]
                }
                ss_data[Pos1<Pos2 &
                          pos_diag == guide$pos_guide[k] &
                          pos_perp %in% helper[!is.na(p_strand),pos],
                        sheet_antipar_nr:=helper[pos==pos_perp,min(closeby_strands)],pos_perp]
                ss_data[Pos1<Pos2 & sheet_antipar_nr == min(closeby_strands),
                        sheet_antipar_p := sumlog(beta_antipar_p) $p]
              } else {
                max_strand = max(c(ss_data$sheet_antipar_nr,0),na.rm=T)
                ss_data[Pos1<Pos2 &
                          pos_diag == guide$pos_guide[k] &
                          pos_perp %in% helper[!is.na(p_strand),pos],
                        sheet_antipar_nr:=helper[pos==pos_perp,strand + max_strand],pos_perp]
                
                ss_data[Pos1<Pos2 &
                          pos_diag == guide$pos_guide[k] &
                          pos_perp %in% helper[!is.na(p_strand),pos],
                        sheet_antipar_p:=helper[pos==pos_perp,p_strand],pos_perp]
              }
            }
          }
        }
      }
      return(ss_data)
    }
    ##################################################
    
    ss_data = calculate_seed_pval(ss_data,restricted_pairing,input_ss,seed_size)
    
    
    #######################################
    ######### plot initial results ######## 
    #######################################
    
    ###raw data
    PLOT = ggplot() +
      geom_raster(data=ss_data,aes(Pos1,Pos2,fill=input)) +
      geom_raster(data=ss_data,aes(Pos2,Pos1,fill=input)) +
      scale_fill_gradient2(midpoint=0,low="tomato3",high="steelblue3",na.value = "white") +
      
      
      geom_abline(linetype=3,alpha=0.33,slope = -1) +
      scale_x_continuous(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
      scale_y_reverse(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
      labs(fill="input")
    #alpha helices
    if (nrow(input_ss[SS=="H"]) > 0) {
      PLOT = PLOT + geom_rect(data=input_ss[SS=="H",.(pos_min = min(Pos),pos_max = max(Pos)),rleidx],
                              inherit.aes = F,
                              aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="darkgreen",alpha=0)
    }
    #beta strands
    if (nrow(input_ss[SS=="E"]) > 0) {
      PLOT = PLOT + geom_rect(data=input_ss[SS=="E",.(pos_min = min(Pos),pos_max = max(Pos)),rleidx],
                              inherit.aes = F,
                              aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="red",alpha=0)
      #potential beta sheets
      if (nrow(input_ss[SS=="E"]) > 1) {
        helper1 = expand.grid(unique(input_ss[SS == "E",rleidx]),unique(input_ss[SS == "E",rleidx]))
        potential_beta_sheets = data.table(s1=helper1[,1],s2=helper1[,2])
        potential_beta_sheets[,':=' (xmin = input_ss[rleidx==s1,min(Pos)-0.5],xmax = input_ss[rleidx==s1,max(Pos)+0.5]),s1]
        potential_beta_sheets[,':=' (ymin = input_ss[rleidx==s2,min(Pos)-0.5],ymax = input_ss[rleidx==s2,max(Pos)+0.5]),s2]
        PLOT = PLOT + geom_rect(data=potential_beta_sheets[s1!=s2],inherit.aes=F,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),color="orange",linetype=2,alpha=0)
      }
    }
    ggsave(plot = PLOT, file=paste0(dataset_dir,"results/secondary_structure/",prefix,eval_cols[eval_cols_idx],"_betasheet_input_rawscores.pdf"),width=6,height=5)
    
    ###sheet propensity after smoothing
    PLOT = ggplot() +
      geom_raster(data=ss_data[Pos1 < Pos2],aes(Pos1,Pos2,fill=log10(beta_par_p))) +
      geom_raster(data=ss_data[Pos1 < Pos2],aes(Pos2,Pos1,fill=log10(beta_antipar_p))) +
      scale_fill_distiller(na.value="white") +
      geom_abline(slope=-1,linetype=3) +
      geom_abline(linetype=3,alpha=0.33,slope = -1) +
      scale_x_continuous(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
      scale_y_reverse(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
      labs(fill="log10(p)")
    #alpha helices
    if (nrow(input_ss[SS=="H"]) > 0) {
      PLOT = PLOT + geom_rect(data=input_ss[SS=="H",.(pos_min = min(Pos),pos_max = max(Pos)),rleidx],
                              inherit.aes = F,
                              aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="darkgreen",alpha=0)
    }
    #beta strands
    if (nrow(input_ss[SS=="E"]) > 0) {
      PLOT = PLOT + geom_rect(data=input_ss[SS=="E",.(pos_min = min(Pos),pos_max = max(Pos)),rleidx],
                              inherit.aes = F,
                              aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="red",alpha=0)
      #potential beta sheets
      if (nrow(input_ss[SS=="E"]) > 1) {
        helper1 = expand.grid(unique(input_ss[SS == "E",rleidx]),unique(input_ss[SS == "E",rleidx]))
        potential_beta_sheets = data.table(s1=helper1[,1],s2=helper1[,2])
        potential_beta_sheets[,':=' (xmin = input_ss[rleidx==s1,min(Pos)-0.5],xmax = input_ss[rleidx==s1,max(Pos)+0.5]),s1]
        potential_beta_sheets[,':=' (ymin = input_ss[rleidx==s2,min(Pos)-0.5],ymax = input_ss[rleidx==s2,max(Pos)+0.5]),s2]
        PLOT = PLOT + geom_rect(data=potential_beta_sheets[s1!=s2],inherit.aes=F,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),color="orange",linetype=2,alpha=0)
      }
    }
    ggsave(plot = PLOT,file=paste0(dataset_dir,"results/secondary_structure/",prefix,eval_cols[eval_cols_idx],"_betasheet_input_sheet_smooth.pdf"),width=6,height=5)
    
    ### propensities of beta sheet seeds
    PLOT = ggplot() +
      geom_raster(data=ss_data[Pos1<Pos2],aes(Pos1,Pos2,fill=log10(sheet_par_p))) +
      geom_raster(data=ss_data[Pos1<Pos2],aes(Pos2,Pos1,fill=log10(sheet_antipar_p))) +
      scale_fill_distiller(na.value="white") +
      geom_abline(linetype=3,alpha=0.33,slope = -1) +
      scale_x_continuous(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
      scale_y_reverse(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
      labs(fill="log10(p)")
    #alpha helices
    if (nrow(input_ss[SS=="H"]) > 0) {
      PLOT = PLOT + geom_rect(data=input_ss[SS=="H",.(pos_min = min(Pos),pos_max = max(Pos)),rleidx],
                              inherit.aes = F,
                              aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="darkgreen",alpha=0)
    }
    #beta strands
    if (nrow(input_ss[SS=="E"]) > 0) {
      PLOT = PLOT + geom_rect(data=input_ss[SS=="E",.(pos_min = min(Pos),pos_max = max(Pos)),rleidx],
                              inherit.aes = F,
                              aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="red",alpha=0)
      #potential beta sheets
      if (nrow(input_ss[SS=="E"]) > 1) {
        helper1 = expand.grid(unique(input_ss[SS == "E",rleidx]),unique(input_ss[SS == "E",rleidx]))
        potential_beta_sheets = data.table(s1=helper1[,1],s2=helper1[,2])
        potential_beta_sheets[,':=' (xmin = input_ss[rleidx==s1,min(Pos)-0.5],xmax = input_ss[rleidx==s1,max(Pos)+0.5]),s1]
        potential_beta_sheets[,':=' (ymin = input_ss[rleidx==s2,min(Pos)-0.5],ymax = input_ss[rleidx==s2,max(Pos)+0.5]),s2]
        PLOT = PLOT + geom_rect(data=potential_beta_sheets[s1!=s2],inherit.aes=F,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),color="orange",linetype=2,alpha=0)
      }
    }
    ggsave(plot = PLOT,file=paste0(dataset_dir,"results/secondary_structure/",prefix,eval_cols[eval_cols_idx],"_betasheet_input_sheet_seeds.pdf"),width=6,height=5)
    
    
    
    
    
    ###########################################
    ######## decide beta sheet pairing ########
    ###########################################
    
    
    ##########
    ## function to plot intermediate states of the process
    plot_intermediate_function  = function(plot_int_idx,ss_data,beta_strands,interaction_data,interaction_helper) {
      
      # browser()
      
      PLOT = ggplot() +
        geom_raster(data=ss_data[Pos1<=Pos2],aes(Pos1,Pos2,fill=log10(beta_par_p))) +
        geom_raster(data=ss_data[Pos1<Pos2],aes(Pos2,Pos1,fill=log10(beta_antipar_p))) +
        scale_fill_distiller(na.value="white") +
        geom_abline(linetype=3,alpha=0.33,slope = -1) +
        scale_x_continuous(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
        scale_y_reverse(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
        labs(fill="log10(p)")
      #ss from input_ss
      #alpha helices
      if (nrow(input_ss[SS=="H"]) > 0) {
        PLOT = PLOT + geom_rect(data=input_ss[SS=="H",.(pos_min = min(Pos),pos_max = max(Pos)),rleidx],
                                inherit.aes = F,
                                aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="red",alpha=0)
      }
      #beta strands
      if (nrow(input_ss[SS=="E"]) > 0) {
        PLOT = PLOT + geom_rect(data=input_ss[SS=="E",.(pos_min = min(Pos),pos_max = max(Pos)),rleidx],
                                inherit.aes = F,
                                aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="darkgreen",alpha=0)
      }
      #beta sheet interactions
      if (nrow(interaction_data)>0) {
        PLOT = PLOT + geom_rect(data=interaction_data[type=="par"],aes(xmin=pos1_min-.5,xmax=pos1_max+.5,ymin=pos2_min-.5,ymax=pos2_max+.5),
                                color="orange",alpha=0)
        PLOT = PLOT + geom_rect(data=interaction_data[type=="anti-par"],aes(xmin=pos2_min-.5,xmax=pos2_max+.5,ymin=pos1_min-.5,ymax=pos1_max+.5),
                                color="orange",alpha=0)
      }
      #current interaction that causes reset
      if (nrow(interaction_helper)>0) {
        if (interaction_helper$type == "par") {
          PLOT = PLOT + geom_rect(data=interaction_helper,aes(xmin=pos1_min-.5,xmax=pos1_max+.5,ymin=pos2_min-.5,ymax=pos2_max+.5),
                                  color="black",alpha=0)
        } else {
          PLOT = PLOT + geom_rect(data=interaction_helper,aes(xmin=pos2_min-.5,xmax=pos2_max+.5,ymin=pos1_min-.5,ymax=pos1_max+.5),
                                  color="black",alpha=0)
        }
      }
      
      ggsave(plot = PLOT,file=paste0(dataset_dir,"results/secondary_structure/temp_plots/",prefix,eval_cols[eval_cols_idx],"secondary_structure_map_idx",plot_int_idx,".pdf"),width=6,height=5)
      plot_int_idx = plot_int_idx+1
      return(plot_int_idx)
    }
    ##########
    
    ##########
    #function to calculate overlaps of potential beta sheets with known strands and blocked space
    compute_strand_overlap = function(potential_sheets,beta_strands,blocked_space) {
      #get strand overlaps
      setkey(potential_sheets,pos1_min,pos1_max)
      setkey(beta_strands,pos_min,pos_max)
      
      #record lowest(meaning lowest position!) strand number overlap
      potential_sheets[,strand1 := foverlaps(potential_sheets,beta_strands,type="any", mult="first")[sheet_nr == unlist(.SD),strand_nr],
                       sheet_nr,.SDcols = "sheet_nr"]
      #record number of overlaps
      potential_sheets[,overlaps1 := foverlaps(potential_sheets,beta_strands,type="any")[sheet_nr == unlist(.SD),.N],
                       sheet_nr,.SDcols = "sheet_nr"]
      #if overlapping with alpha-helix, delete
      potential_sheets[sheet_nr %in% foverlaps(potential_sheets,blocked_space,type="any", mult="first")[!is.na(helix_nr),sheet_nr],strand1 := NA]
      #same for position 2
      setkey(potential_sheets,pos2_min,pos2_max)
      potential_sheets[,strand2 := foverlaps(potential_sheets,beta_strands,type="any", mult="first")[sheet_nr == unlist(.SD),strand_nr],
                       sheet_nr,.SDcols = "sheet_nr"]
      potential_sheets[,overlaps2 := foverlaps(potential_sheets,beta_strands,type="any")[sheet_nr == unlist(.SD),.N],
                       sheet_nr,.SDcols = "sheet_nr"]
      potential_sheets[sheet_nr %in% foverlaps(potential_sheets,blocked_space,type="any", mult="first")[!is.na(helix_nr),sheet_nr],strand2 := NA]
      setkey(potential_sheets,pval)
      return(potential_sheets)
    }
    ##########
    
    ##########
    #function to calculate overlaps of potential beta sheets with known strands and blocked space
    compute_strand_overlap_v2 = function(ss_data,beta_strands,blocked_space) {
      
      #calculate seed pvalues
      if (restricted_pairing == F) {
        ss_data = calculate_seed_pval(ss_data = ss_data[,.SD,,.SDcols = setdiff(1:ncol(ss_data),grep("sheet_",names(ss_data)))],
                                      restricted_pairing = F,input_ss = c(),
                                      blocked_space = blocked_space)
      }
      # if (debug_this) {browser()}
      
      direction = c("antipar","par")
      helper = list()
      for (i in 1:2) {
        if (i == 1) {
          potential_sheets = unique(ss_data[!is.na(sheet_antipar_nr),
                                            .(pval = sheet_antipar_p,
                                              pos1_min = min(Pos1),pos1_max = max(Pos1),pos2_min = min(Pos2),pos2_max= max(Pos2),
                                              type = "anti-par"),sheet_antipar_nr])[order(pval)]
          names(potential_sheets)[grep(names(potential_sheets),pattern="sheet_antipar_nr")] = "sheet_nr"
        } else {
          potential_sheets = unique(ss_data[!is.na(sheet_par_nr),
                                            .(pval = sheet_par_p,
                                              pos1_min = min(Pos1),pos1_max = max(Pos1),pos2_min = min(Pos2),pos2_max = max(Pos2),
                                              type = "par"),sheet_par_nr])[order(pval)]
          names(potential_sheets)[grep(names(potential_sheets),pattern="sheet_par_nr")] = "sheet_nr"
        }
        
        #get strand overlaps
        setkey(potential_sheets,pos1_min,pos1_max)
        setkey(beta_strands,pos_min,pos_max)
        
        #record lowest(meaning lowest position!) strand number overlap
        potential_sheets[,strand1 := foverlaps(potential_sheets,beta_strands,type="any", mult="first")[sheet_nr == unlist(.SD),strand_nr],
                         sheet_nr,.SDcols = "sheet_nr"]
        #record number of overlaps
        potential_sheets[,overlaps1 := foverlaps(potential_sheets,beta_strands,type="any")[sheet_nr == unlist(.SD),.N],
                         sheet_nr,.SDcols = "sheet_nr"]
        #if overlapping with alpha-helix, delete
        potential_sheets[sheet_nr %in% foverlaps(potential_sheets,blocked_space,type="any", mult="first")[!is.na(helix_nr),sheet_nr],strand1 := NA]
        #same for position 2
        setkey(potential_sheets,pos2_min,pos2_max)
        potential_sheets[,strand2 := foverlaps(potential_sheets,beta_strands,type="any", mult="first")[sheet_nr == unlist(.SD),strand_nr],
                         sheet_nr,.SDcols = "sheet_nr"]
        potential_sheets[,overlaps2 := foverlaps(potential_sheets,beta_strands,type="any")[sheet_nr == unlist(.SD),.N],
                         sheet_nr,.SDcols = "sheet_nr"]
        potential_sheets[sheet_nr %in% foverlaps(potential_sheets,blocked_space,type="any", mult="first")[!is.na(helix_nr),sheet_nr],strand2 := NA]
        setkey(potential_sheets,pval)
        
        if (i == 1) {
          helper[[1]] = copy(potential_sheets)
        } else {
          helper[[2]] = copy(potential_sheets)
        }
      }
      return(helper)
    }
    ##########
    
    ##########
    decide_beta_pairing = function(ss_data,input_ss,integration_threshold = 10^-3,trim_beta_strand = TRUE,plot_intermediate = FALSE,debug_this=FALSE) {
      
      ### record known beta strands
      beta_strands = input_ss[SS == "E",.(pos_min = min(Pos),pos_max = max(Pos),interactions = 0,pval = as.numeric(1)),rleidx]
      beta_strands[,strand_nr := 1:.N]
      beta_strands[,rleidx := NULL]
      setkey(beta_strands,pos_min,pos_max)
      
      ### define space blocked for new beta strands and their interactions
      #first known alpha helices, later also space between split beta strands
      if (nrow(input_ss[SS=="H"]) > 0) {
        blocked_space = input_ss[SS=="H",.(pos_min = min(Pos),pos_max = max(Pos),helix_nr = rleidx),rleidx]
        #also block up to three positions up and downstream of helices
        for (b in nrow(blocked_space)) {
          up = T
          down = T
          for (idx in 1:3) {
            if (down ==T & input_ss[Pos == blocked_space[b,pos_min-1],SS == "C"]) {
              blocked_space[b,pos_min := pos_min - 1]
            } else {
              down == F
            }
            if (up ==T & input_ss[Pos == blocked_space[b,pos_max+1],SS == "C"]) {
              blocked_space[b,pos_max := pos_max + 1]
            } else {
              up == F
            }
          }
        }
        blocked_space[,rleidx := NULL]
      } else {
        blocked_space = data.table(helix_nr = 0, pos_min = 0, pos_max = 0)
      }
      setkey(blocked_space,pos_min,pos_max)
      
      ### find potential beta sheet pairing in anti-parallel direction
      helper = compute_strand_overlap_v2(ss_data,beta_strands,blocked_space)
      lowest_p_antipar = helper[[1]]
      lowest_p_par = helper[[2]]
      
      ### initialize interaction_data
      interaction_data = data.table(sheet_nr = numeric(),pval = numeric(),pos1_min = numeric(),pos1_max = numeric(),
                                    pos2_min = numeric(),pos2_max = numeric(),type = character(),
                                    strand1 = integer(),overlaps1 = numeric(),strand2 = integer(),overlaps2 = numeric())
      
      if (debug_this) {
        print("initial setup")
        print(beta_strands)
        print(blocked_space)
        print(lowest_p_par)
        print(lowest_p_antipar)
        browser()
      }
      
      ### plotting intemediate results?
      if (plot_intermediate == TRUE) {plot_int_idx = 1} else {plot_int_idx = 0}
      
      ######### run sheet detection
      while (min(c(lowest_p_antipar[1]$pval,lowest_p_par[1]$pval),na.rm=T) < integration_threshold) {
        
        reset = FALSE
        add_interaction = TRUE
        
        #choose anti-par or par beta sheet interaction
        if (which.min(c(lowest_p_antipar[1]$pval,lowest_p_par[1]$pval)) == 1) { #anti-par
          interaction_helper = lowest_p_antipar[1] } else {
            interaction_helper = lowest_p_par[1] } #par
        
        if (debug_this) {
          print("new interaction_helper")
          print(interaction_helper)
          browser()
        }
        
        if (is.na(interaction_helper$strand1) | is.na(interaction_helper$strand2)) { #if interaction has at least one strand not detected previously
          
          if (is.na(interaction_helper$strand1) & is.na(interaction_helper$strand2)) {
            add_interaction = FALSE
          } else if (is.na(interaction_helper$strand1)) {
            #check whether the unknown strand overlaps with an alpha helix
            if (foverlaps(interaction_helper,blocked_space,by.x=c("pos1_min","pos1_max"))[,sum(!is.na(helix_nr))>0]) {
              #does overlap with existing alpha helix -> delete interaction
              add_interaction = FALSE
            } else if (beta_strands[strand_nr == interaction_helper$strand2,interactions == 2]) { #does other strand already have 2 interactions?
              add_interaction = FALSE
            } else {
              
              if (plot_intermediate) {
                plot_int_idx = plot_intermediate_function(plot_int_idx = plot_int_idx,ss_data = ss_data,beta_strands = beta_strands,
                                                          interaction_data = interaction_data,interaction_helper =  interaction_helper)
              }
              
              ## modify beta strands
              beta_strands = rbind(beta_strands,data.table(strand_nr = as.integer(nrow(beta_strands)+1),
                                                           pval = sumlog(ss_data[Pos1==Pos2 & between(Pos1,interaction_helper$pos1_min,interaction_helper$pos1_max),beta_par_p])$p,
                                                           pos_min = interaction_helper$pos1_min,
                                                           pos_max = interaction_helper$pos1_max,
                                                           interactions = 0))
              #reset strand_nrs
              setorder(beta_strands,pos_min); beta_strands[,strand_nr := 1:.N]
              
              if (debug_this) {
                print("new beta strand given strand 1 NA")
                print(beta_strands)
                browser()
              }
              
              reset = TRUE
            }
          } else if (is.na(interaction_helper$strand2)) {
            if (foverlaps(interaction_helper,blocked_space,by.x=c("pos2_min","pos2_max"))[,sum(!is.na(helix_nr))>0]) {
              #does overlap with existing alpha helix -> delete interaction
              add_interaction = FALSE
            } else if (beta_strands[strand_nr == interaction_helper$strand1,interactions == 2]) { #does other strand already have 2 interactions?
              add_interaction = FALSE
            } else {
              
              if (plot_intermediate) {
                plot_int_idx = plot_intermediate_function(plot_int_idx = plot_int_idx,ss_data = ss_data,beta_strands = beta_strands,
                                                          interaction_data = interaction_data,interaction_helper =  interaction_helper)
              }
              
              ## modify beta strands
              beta_strands = rbind(beta_strands,data.table(strand_nr = as.integer(nrow(beta_strands)+1),
                                                           pval = sumlog(ss_data[Pos1==Pos2 & between(Pos1,interaction_helper$pos2_min,interaction_helper$pos2_max),beta_par_p])$p,
                                                           pos_min = interaction_helper$pos2_min,
                                                           pos_max = interaction_helper$pos2_max,
                                                           interactions = 0))
              #reset strand_nrs
              setorder(beta_strands,pos_min); beta_strands[,strand_nr := 1:.N]
              
              if (debug_this) {
                print("new beta strand given strand 2 NA")
                print(beta_strands)
                browser()
              }
              
              reset = TRUE
            }
          }
        } else if (interaction_helper$strand1 == interaction_helper$strand2) { #if it is an intra-strand interaction, split strand
          old_strand = beta_strands[strand_nr == interaction_helper$strand1]
          
          if (interaction_helper$type == "anti-par") {
            if (plot_intermediate) {
              plot_int_idx = plot_intermediate_function(plot_int_idx = plot_int_idx,ss_data = ss_data,beta_strands = beta_strands,
                                                        interaction_data = interaction_data,interaction_helper =  interaction_helper)
            }
            #define breakpoint
            strand_break = interaction_helper[,(pos1_max + pos2_min)/2]
            
            ##### check whether proposed strandbreak is not already covered by an interaction
            helper = data.table(pos_min = ceiling(strand_break-1),pos_max = floor(strand_break+1))
            setkey(helper,pos_min,pos_max)
            setkey(interaction_data,pos1_min,pos1_max)
            pos1 = foverlaps(interaction_data,helper)
            setkey(interaction_data,pos2_min,pos2_max)
            pos2 = foverlaps(interaction_data,helper)
            if (sum(!is.na(pos1$pos_min)) > 0 | sum(!is.na(pos2$pos_min)) > 0) { #proposed strandbreak overlaps with already established interaction
              add_interaction = FALSE
            } else {
              new_strands = rbind(old_strand,old_strand)
              new_strands[1,pos_max := ceiling(strand_break-2)]
              #update pval
              new_strands[1,pval := sumlog(ss_data[Pos1==Pos2 & between(Pos1,pos_min,pos_max),beta_par_p])$p]
              new_strands[2,strand_nr := strand_nr+1]
              new_strands[2,pos_min := floor(strand_break+2)]
              new_strands[2,pval := sumlog(ss_data[Pos1==Pos2 & between(Pos1,pos_min,pos_max),beta_par_p])$p]
              
              #update all strands
              beta_strands[strand_nr > old_strand$strand_nr,strand_nr := strand_nr + 1L]
              
              beta_strands = rbind(beta_strands[strand_nr != old_strand$strand_nr],new_strands)
              #reset strand_nrs
              setorder(beta_strands,pos_min); beta_strands[,strand_nr := as.integer(1:.N)]
              setkey(beta_strands,pos_min,pos_max)
              
              ## define blocked space where strands were split ##
              blocked_space = rbind(blocked_space,data.table(helix_nr = 666,pos_min = ceiling(strand_break-1), pos_max = floor(strand_break+1)))
              setkey(blocked_space,pos_min,pos_max)
              
              if (debug_this) {
                print("strand split")
                print(beta_strands)
                print(blocked_space)
                browser()
              }
              
              reset = TRUE
            }
          } else {
            if (plot_intermediate) {
              plot_int_idx = plot_int_idx = plot_intermediate_function(plot_int_idx = plot_int_idx,ss_data = ss_data,beta_strands = beta_strands,
                                                                       interaction_data = interaction_data,interaction_helper =  interaction_helper)
            }
            print("help, splitting strands if interaction == par not defined yet")
            if (debug_this) {browser()}
          }
        }
        
        if (debug_this) {
          print(paste0("reset: ",as.character(reset)));
          print(paste0("add interaction: ",as.character(add_interaction)))
          browser()
        }
        
        
        if (reset == TRUE) {
          ### reset the whole process, easier than tinkering with everything upstream
          beta_strands[,interactions := 0]
          setkey(beta_strands,pos_min,pos_max)
          
          ### recompute potential beta sheet pairing
          helper = compute_strand_overlap_v2(ss_data,beta_strands,blocked_space)
          lowest_p_antipar = helper[[1]]
          lowest_p_par = helper[[2]]
          
          #re-set interaction_data
          interaction_data = data.table(sheet_nr = numeric(),pval = numeric(),pos1_min = numeric(),pos1_max = numeric(),
                                        pos2_min = numeric(),pos2_max = numeric(),type = character(),
                                        strand1 = numeric(),overlaps1 = numeric(),strand2 = numeric(),overlaps2 = numeric())
          
          
        } else if (add_interaction == TRUE) { #if interaction falls within two separate strands, extend strands if necessary and shorten other strand
          
          if (beta_strands[strand_nr == interaction_helper$strand1,interactions] < 2 &
              beta_strands[strand_nr == interaction_helper$strand2,interactions] < 2 &
              sum(beta_strands$interactions) < (nrow(beta_strands) + max(0,(nrow(beta_strands)-2)))) { #if both strands have fewer than 2 interactions, add the next one
            #else: just delete this interaction and all interactions for same strands
            
            #strand1
            beta_strands[strand_nr == interaction_helper$strand1,pos_min := min(pos_min,interaction_helper$pos1_min)]
            beta_strands[strand_nr == interaction_helper$strand1,
                         pos_max := min(max(pos_max,interaction_helper$pos1_max),interaction_helper$pos2_min-3)]
            
            #strand2
            beta_strands[strand_nr == interaction_helper$strand2,
                         pos_min := max(min(pos_min,interaction_helper$pos2_min),interaction_helper$pos1_max+3)]
            beta_strands[strand_nr == interaction_helper$strand2,pos_max := max(pos_max,interaction_helper$pos2_max)]
            
            #recompute association of potential sheet interacitons with strands (after updating strand positions)
            helper = compute_strand_overlap_v2(ss_data,beta_strands,blocked_space)
            lowest_p_antipar = helper[[1]]
            lowest_p_par = helper[[2]]
            # lowest_p_antipar = compute_strand_overlap(lowest_p_antipar,beta_strands,blocked_space)
            # lowest_p_par = compute_strand_overlap(lowest_p_par,beta_strands,blocked_space)
            
            beta_strands[strand_nr == interaction_helper$strand1,interactions := interactions + 1]
            beta_strands[strand_nr == interaction_helper$strand2,interactions := interactions + 1]
            
            if (nrow(interaction_data) > 0) {
              interaction_data = rbind(interaction_data,interaction_helper)
            } else { interaction_data = interaction_helper }
            
            #delete all sheet interactions that would fall within same strand pair
            for (i in 1:nrow(interaction_data)) {
              lowest_p_antipar = lowest_p_antipar[!(strand1 == interaction_data[i]$strand1 & strand2 == interaction_data[i]$strand2) | is.na(strand1) | is.na(strand2)]
              lowest_p_par = lowest_p_par[!(strand1 == interaction_data[i]$strand1 & strand2 == interaction_data[i]$strand2) | is.na(strand1) | is.na(strand2)]  
            }
            
          } else {
            add_interaction = FALSE
            if (debug_this) {
              print("add_interaction aborted, too many interactions on strands already")
              browser()}
          }
        }
        
        #delete specific interaction if it doesn't fit in
        if (add_interaction == FALSE & interaction_helper$type == "anti-par") {
          lowest_p_antipar = lowest_p_antipar[2:.N]
        } else if (add_interaction == FALSE & interaction_helper$type == "par"){
          lowest_p_par = lowest_p_par[2:.N]
        }
        
        if (debug_this) {
          print("finished iteration")
          print(beta_strands)
          print(blocked_space)
          if (exists("interaction_data")) {print(interaction_data)}
          browser()
        }
        
      }
      
      ### print results
      print(beta_strands)
      print(blocked_space)
      if (exists("interaction_data")) {print(interaction_data)}
      
      ### transfer results
      ss_data[,strand_final_nr := as.numeric(NA)]
      ss_data[,strand_final_p := as.numeric(NA)]
      
      if (exists("interaction_data")) {
        if (trim_beta_strand) { #if TRUE, define only as beta strand what pairs as a betasheet
          strands = unique(c(interaction_data$strand1,interaction_data$strand2))
          for (i in strands) {
            helper = interaction_data[strand1 == i,.(pos_min = pos1_min,pos_max = pos1_max)]
            helper = rbind(helper,interaction_data[strand2 == i,.(pos_min = pos2_min,pos_max = pos2_max)])
            ss_data[Pos1 == Pos2 & between(Pos1,helper[,min(pos_min)],helper[,max(pos_max)]),
                    ':=' (strand_final_nr = i,strand_final_p = sumlog(beta_par_p)$p)]
          }
        } else { #if FALSE, keep beta strand overlaps without betasheet pairing
          for (i in 1:nrow(beta_strands)) {
            ss_data[Pos1 == Pos2 & between(Pos1,beta_strands[i,pos_min],beta_strands[i,pos_max]),
                    ':=' (strand_final_nr = beta_strands[i,strand_nr],strand_final_p = beta_strands[i,pval])]
          }
        }
      }
      if (plot_intermediate) {
        plot_int_idx = plot_intermediate_function(plot_int_idx = 999,ss_data = ss_data,beta_strands = beta_strands,
                                                  interaction_data = interaction_data,interaction_helper =  data.table(bla = integer()))
      }
      
      if (debug_this) {browser()}
      
      ### return results
      if (exists("interaction_data")) {
        return_list = list(ss_data,interaction_data)
      } else {
        return_list = list(ss_data,c())
      }
      return(return_list)
    }
    #########################################################################################################
    
    
    A = decide_beta_pairing(ss_data,
                            input_ss,
                            trim_beta_strand = !restricted_pairing,
                            plot_intermediate = TRUE, 
                            debug_this = debug_this,
                            integration_threshold=0.05)
    ss_data = A[[1]]
    beta_sheet_pairing = A[[2]]
    
    ss_elements = data.table(position = pos_range[1]:pos_range[2],ss = "C")
    ss_elements[position %in% ss_data[Pos1 == Pos2 & !is.na(strand_final_nr),Pos1],ss := "E"]
    ss_elements[position %in% input_ss[SS == "H",Pos],ss := "H"]
    
    ######################
    #### save results ####
    ######################
    save(file=paste0(dataset_dir,"processed_data/",prefix,"secondary_structure_elements_",eval_cols[eval_cols_idx],".RData"),
         list = c("ss_data","beta_sheet_pairing","ss_elements"))
    
    
    
    ########################################
    #### plot final predictions results ####
    ########################################
    
    #get known secondary structure elements to compare to
    if (!is.null(known_ss_file)) {
      known_ss = fread(known_ss_file)
      names(known_ss)[2] = "SS"
      known_ss[,rleidx := rleid(SS)]
    } else {
      known_ss = copy(input_ss)
    }
    
    ### get known beta sheet interactions
    if (!is.null(known_bsi_file)) {
      known_beta_sheets = fread(known_bsi_file)
      known_beta_sheets[pos_o > pos_hn,c("pos_o","pos_hn") := .(pos_hn,pos_o)] # for plotting par/anti-par hbonds consistently on one side diagonal
      
      known_beta_sheets_blocks = data.table(pos1_min = as.numeric(),pos1_max = as.numeric(),
                                            pos2_min = as.numeric(),pos2_max = as.numeric(), 
                                            type = as.character(),strand1 = as.numeric(),strand2 = as.numeric())
      if (nrow(known_beta_sheets[type=="anti_par"]) > 0) {
        known_beta_sheets_blocks = rbind(known_beta_sheets_blocks,
                                         known_beta_sheets[type=="anti_par",.(pos1_min = min(pos_hn),pos1_max = max(pos_hn),
                                                                              pos2_min = min(pos_o),pos2_max = max(pos_o), type = "anti_par"),.(strand1,strand2)])
      }
      
      if (nrow(known_beta_sheets[type=="par"]) > 0) {
        par_strands_unique = unique(known_beta_sheets[type=="par",.(strand1,strand2)])
        for (i in 1:nrow(par_strands_unique)) {
          helper = known_beta_sheets[type == "par" & strand1==par_strands_unique[i]$strand1 & strand2==par_strands_unique[i]$strand2]
          known_beta_sheets_blocks = rbind(known_beta_sheets_blocks,
                                           data.table(strand1 = par_strands_unique[i]$strand1,
                                                      strand2 = par_strands_unique[i]$strand2,
                                                      pos1_min = helper[,.N,pos_hn][1,ifelse(N==1,pos_hn,pos_hn-1)],
                                                      pos1_max = helper[,.N,pos_hn][.N,ifelse(N==1,pos_hn,pos_hn+1)],
                                                      pos2_min = as.integer(helper[,.N,pos_o][1,ifelse(N==1,pos_o,pos_o-1)]),
                                                      pos2_max = as.integer(helper[,.N,pos_o][.N,ifelse(N==1,pos_o,pos_o+1)]), type="par"))
        }
        known_beta_sheets_blocks[type == "par" & pos1_min > pos2_min, c("pos1_min","pos1_max","pos2_min","pos2_max") := .(pos2_min,pos2_max,pos1_min,pos1_max)]
      }
    } else {
      known_beta_sheets = data.table(pos_hn = as.numeric(), pos_o = as.numeric(), strand1 = as.numeric(), strand2=as.numeric(), type=as.character())
    }
    
    ##################
    ## plot final data with secondary structures and beta sheet pairing
    PLOT = ggplot() +
      geom_raster(data=ss_data[Pos1<Pos2],aes(Pos1,Pos2,fill=log10(sheet_par_p))) +
      geom_raster(data=ss_data[Pos1==Pos2],aes(Pos1,Pos2,fill=log10(pmin(strand_final_p,na.rm=T)))) +
      geom_raster(data=ss_data[Pos1<Pos2],aes(Pos2,Pos1,fill=log10(sheet_antipar_p))) +
      scale_fill_distiller(na.value="white",limits=c(ss_data[,floor(log10(min(sheet_par_p,strand_final_p,sheet_antipar_p,na.rm=T)))],0)) +
      geom_abline(linetype=3,alpha=0.33,slope = -1) +
      scale_x_continuous(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
      scale_y_reverse(breaks = seq(5,pos_range[2],5),expand = c(0.01,0)) +
      labs(fill="log10(p)")
    
    #ss from DMS
    #alpha helices
    if (nrow(input_ss[SS=="H"]) > 0) {
      PLOT = PLOT + geom_rect(data=input_ss[SS=="H",.(pos_min = min(Pos),pos_max = max(Pos)),rleidx],
                              inherit.aes = F,
                              aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="darkgreen",alpha=0)
    }
    #beta strands (updated from beta sheet pairing)
    if (nrow(input_ss[SS=="E"]) > 0) {
      PLOT = PLOT + geom_rect(data=ss_data[Pos1==Pos2 & !is.na(strand_final_nr),.(pos_min = min(Pos1),pos_max = max(Pos1)),strand_final_nr],
                              inherit.aes = F,
                              aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5),color="red",alpha=0)
    }
    
    #beta sheet interactions
    if (nrow(beta_sheet_pairing[type=="par"]) > 0) {
      PLOT = PLOT + geom_rect(data=beta_sheet_pairing[type=="par"],aes(xmin=pos1_min-.5,xmax=pos1_max+.5,ymin=pos2_min-.5,ymax=pos2_max+.5),
                              color="orange",alpha=0) +
        geom_line(data=beta_sheet_pairing[type=="par",.(x=c(pos1_min-.5,pos1_max+.5),y=c(pos2_min-.5,pos2_max+.5)),sheet_nr],aes(x,y,group=sheet_nr),
                  color="orange",alpha=1,linetype=2)
    }
    if (nrow(beta_sheet_pairing[type=="anti-par"]) > 0) {
      PLOT = PLOT + geom_rect(data=beta_sheet_pairing[type=="anti-par"],aes(xmin=pos2_min-.5,xmax=pos2_max+.5,ymin=pos1_min-.5,ymax=pos1_max+.5),
                              color="orange",alpha=0) +
        geom_line(data=beta_sheet_pairing[type=="anti-par",.(x=c(pos2_min-.5,pos2_max+.5),y=c(pos1_max+.5,pos1_min-.5)),sheet_nr],aes(x,y,group=sheet_nr),
                  color="orange",alpha=1,linetype=2)
    }
    
    #ss from pdb file
    if (exists("known_beta_sheets_blocks")) {
      PLOT = PLOT + geom_rect(data=known_ss[,.(ss=unique(SS),pos_min = min(Pos),pos_max = max(Pos)),rleid(SS)][ss %in% c("H","E")],
                              inherit.aes=F,
                              aes(xmin=pos_min-.5,xmax=pos_max+.5,ymin=pos_min-.5,ymax=pos_max+.5,linetype=ss),color="black",alpha=0) +
        geom_point(data=known_beta_sheets[type=="par"],aes(pos_o,pos_hn),shape = 1,color="black") +
        geom_point(data=known_beta_sheets[type=="anti_par"],aes(pos_hn,pos_o),shape = 1,color="black") +
        geom_rect(data = known_beta_sheets_blocks,inherit.aes = F,aes(xmin = pos1_min-0.5,xmax = pos1_max + 0.5,ymin = pos2_min-0.5,ymax = pos2_max+0.5),
                  linetype=3,color="black",alpha=0)
    }
    #save
    ggsave(plot = PLOT, file=paste0(dataset_dir,"results/secondary_structure/",prefix,eval_cols[eval_cols_idx],"_betasheet_final.pdf"),width=6,height=5)
    
    
    # browser()
    #precision-recall of predictions
    if (exists("known_beta_sheets_blocks") & nrow(beta_sheet_pairing) > 0) {
      known_beta_sheets_matrix = matrix(0,nrow=pos_range[2],ncol=pos_range[2])
      helper = known_beta_sheets_blocks[pos1_min > pos2_min, c("pos1_min","pos1_max","pos2_min","pos2_max") := .(pos2_min,pos2_max,pos1_min,pos1_max)]
      for (i in 1:nrow(known_beta_sheets_blocks)) {
        if (helper[i,type == "anti_par"]) {
          known_beta_sheets_matrix[cbind(helper[i,seq(pos1_max,pos1_min)],helper[i,seq(pos2_min,pos2_max)])] = 1
        } else {
          known_beta_sheets_matrix[cbind(helper[i,seq(pos2_min,pos2_max)],helper[i,seq(pos1_min,pos1_max)])] = 1
        }
      }
      known_beta_sheets_matrix[1:(pos_range[1]-1),] = 0
      known_beta_sheets_matrix[,1:(pos_range[1]-1)] = 0
      
      beta_sheet_matrix = matrix(0,nrow=pos_range[2],ncol=pos_range[2])
      helper = beta_sheet_pairing[pos1_min > pos2_min, c("pos1_min","pos1_max","pos2_min","pos2_max") := .(pos2_min,pos2_max,pos1_min,pos1_max)]
      for (i in 1:nrow(helper)) {
        if (helper[i,type == "anti-par"]) {
          beta_sheet_matrix[cbind(helper[i,seq(pos1_max,pos1_min)],helper[i,seq(pos2_min,pos2_max)])] = 1
        } else {
          beta_sheet_matrix[cbind(helper[i,seq(pos2_min,pos2_max)],helper[i,seq(pos1_min,pos1_max)])] = 1
        }
      }
      
      prec_recall = data.table(type = "beta sheet", precision = sum(beta_sheet_matrix*known_beta_sheets_matrix)/sum(beta_sheet_matrix),
                               recall = sum(beta_sheet_matrix*known_beta_sheets_matrix)/sum(known_beta_sheets_matrix))
      
      #beta strand
      beta_strand_vec = rep(0,pos_range[2])
      beta_strand_vec[ss_data[Pos1==Pos2 & strand_final_nr,Pos1]] = 1
      known_beta_strand_vec = rep(0,pos_range[2])
      known_beta_strand_vec[known_ss[SS=="E",Pos]] = 1
      known_beta_strand_vec[1:(pos_range[1]-1)] = 0
      prec_recall = rbind(prec_recall,data.table(type = "beta strand", precision = sum(beta_strand_vec*known_beta_strand_vec)/sum(beta_strand_vec),
                                                 recall = sum(beta_strand_vec*known_beta_strand_vec)/sum(known_beta_strand_vec)))
      
      #alpha strand
      alpha_helix_vec = rep(0,pos_range[2])
      alpha_helix_vec[input_ss[SS == "H",Pos]] = 1
      known_alpha_helix_vec = rep(0,pos_range[2])
      known_alpha_helix_vec[known_ss[SS=="H",Pos]] = 1
      known_alpha_helix_vec[1:(pos_range[1]-1)] = 0
      prec_recall = rbind(prec_recall,data.table(type = "alpha helix", precision = sum(alpha_helix_vec*known_alpha_helix_vec)/sum(alpha_helix_vec),
                                                 recall = sum(alpha_helix_vec*known_alpha_helix_vec)/sum(known_alpha_helix_vec)))
      
      prec_recall_melt = melt(prec_recall,id.vars="type")
      PLOT = ggplot(prec_recall_melt,aes(x=type,y=value,fill=variable)) +
        geom_bar(stat="identity",position="dodge") +
        scale_fill_brewer() +
        scale_y_continuous(breaks=seq(0,1,0.1))
      ggsave(plot = PLOT, file=paste0(dataset_dir,"results/secondary_structure/",prefix,eval_cols[eval_cols_idx],"_precision_recall.pdf"),width=6,height=5)
    }
  }
}