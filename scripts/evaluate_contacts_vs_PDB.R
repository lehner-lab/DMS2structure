
#######################################################################################
##### evaluate contacts (top predicted position pairs) against PDB file distances #####
#######################################################################################

evaluate_contacts_vs_PDB = function(contacts,
                                    contactmap,
                                    secondary_structure = NA,
                                    exclude_SS = F,
                                    dataset_dir,
                                    prefix,
                                    N_contactmap = 1,
                                    N_TPrate = 2,
                                    lindist = 5,
                                    dist_cutoff= 8,
                                    modus="cis") {
  
  ### variables 
  # contacts: pairwise interaction score data.table; except for position, wt_aa and #epistatic variant columns, this should only contain the scores that should be evaluated
  # contactmap: contactmap data.table (pairwise positions) from a PDB file; can have multiple distance columns to be evaluated
  # secondary structure: if not NA, data.table with "position" and "ss" identifier column for plotting on contactmaps
  # exclude_SS: if TRUE, position pairs within secondary structure elements are excluded from evaluation
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # N_contactmap: number of top contacts * protein length to plot on contactmaps and use for CDF distance evaluation
  # N_TPrate: number of top contacts * protein length to evaluate for true positive rate (precision) evaluation
  # lindist: only use position pairs with a sequence separation greater than this (only for cis libraries)
  # dist_cutoff: distances cutoff to call a direct contact (in Angstroms)
  # modus: "cis" (single protein) or "trans" (protein-protein interaction)
  
  
  require(data.table)
  require(ggplot2)
  theme_set(theme_minimal(base_size = 9))
  
  eval_cols = setdiff(names(contacts),c("Pos1","Pos2","WT_AA1","WT_AA2","NposE","NnegE"))
  dist_cols = setdiff(names(contactmap),c("Pos1","Pos2","WT_AA1","WT_AA2"))
  
  for (d in seq_along(dist_cols)) {
    data2analyse = merge(contacts,contactmap[,.(Pos1,Pos2,.SD),,.SDcols = dist_cols[d]],by=c("Pos1","Pos2"),all.x = TRUE)
    names(data2analyse)[names(data2analyse) == ".SD"] = "dist"
    
    #tp rate plot
    #ecdf plot
    #contact map for each predictor
    
    #mark regions in data that belong to secondary structure (if a secondary structure data.table is given)
    data2analyse[,ss := FALSE]
    if (!is.na(secondary_structure) & exclude_SS == T & modus == "cis") {
      secondary_structure[,ss_idx := rleid(ss)]
      unique_elements  = secondary_structure[ss %in% c("E","H"),unique(ss_idx)]
      for (i in unique_elements) {
        data2analyse[Pos1 %in% secondary_structure[ss_idx == i,position] & Pos2 %in% secondary_structure[ss_idx == i,position],ss := TRUE]  
      }
    }
    
    
    #define protein length
    protein_length = data2analyse[,diff(range(Pos1))+1] 
    
    # initialize true positive rate (precision) data.table
    TP = data.table(Neval=1:round(N_TPrate*protein_length))
    
    #TP rate for DMS predictions
    for (x in seq_along(eval_cols)) {
      setorderv(data2analyse,eval_cols[x],order = -1,na.last=T)
      if (modus == "cis") {
        TP[,paste0(eval_cols[x]) := data2analyse[Pos1<Pos2-lindist & !is.na(dist) & ss == FALSE &
                                                   !is.na(data2analyse[,unlist(.SD),,.SDcols=eval_cols[x]])
                                                 ][1:round(N_TPrate*protein_length),cumsum(dist<dist_cutoff)/1:round(N_TPrate*protein_length)]]
      } else if (modus == "trans") {
        TP[,paste0(eval_cols[x]) := data2analyse[!is.na(dist) & ss == FALSE &
                                                   !is.na(data2analyse[,unlist(.SD),,.SDcols=eval_cols[x]])
                                                 ][1:round(N_TPrate*protein_length),cumsum(dist<dist_cutoff)/1:round(N_TPrate*protein_length)]]
      }
      
      #CDF of distance of top predicted versus all position pairs
      if (x==1) {
        if (modus == "cis") {
          eCDF = data2analyse[Pos1<Pos2-lindist & !is.na(dist) & ss == FALSE &
                                !is.na(data2analyse[,unlist(.SD),,.SDcols=eval_cols[x]])][1:round(N_contactmap*protein_length),.(dist,id = eval_cols[x])]
        } else if (modus == "trans") {
          eCDF = data2analyse[!is.na(dist) & ss == FALSE &
                                !is.na(data2analyse[,unlist(.SD),,.SDcols=eval_cols[x]])][1:round(N_contactmap*protein_length),.(dist,id = eval_cols[x])]
        }
      } else {
        
        if (modus == "cis") {
          eCDF = rbind(eCDF,data2analyse[Pos1<Pos2-lindist & !is.na(dist) & ss == FALSE &
                                           !is.na(data2analyse[,unlist(.SD),,.SDcols=eval_cols[x]])][1:round(N_contactmap*protein_length),.(dist,id = eval_cols[x])])
        } else if (modus == "trans") {
          eCDF = rbind(eCDF,data2analyse[!is.na(dist) & ss == FALSE &
                                           !is.na(data2analyse[,unlist(.SD),,.SDcols=eval_cols[x]])][1:round(N_contactmap*protein_length),.(dist,id = eval_cols[x])])
        }
      }
    }
    
    if (modus == "cis") {
      eCDF = rbind(eCDF,data2analyse[Pos1<Pos2-lindist,.(dist,id = "xall")])
    } else if (modus == "trans") {
      eCDF = rbind(eCDF,data2analyse[,.(dist,id = "xall")])
    }
    
    measured.vars = eval_cols
    label_cols0 = eval_cols
    
    
    #### write TPrate to file
    write.table(x = round(TP,digits=3), 
                file = paste0(dataset_dir,"results/tertiary_contacts/",prefix,"table_TPrate_lindist",lindist,"_",dist_cols[d],"_",dist_cutoff,"A.text"),
                sep = "\t",row.names = F,quote=F)
    
    
    ##### plot TP rate
    TPmelt = melt(TP,id.vars = c("Neval"),measure.vars = eval_cols)
    
    # calculate random expectation for correct prediction of contacts
    if (modus == "cis") {
      contact_expectation = mean(data2analyse[Pos1<Pos2-lindist,dist<dist_cutoff],na.rm=T)
    } else if (modus == "trans") {
      contact_expectation = mean(data2analyse[,dist<dist_cutoff],na.rm=T)
    }
    
    require(RColorBrewer)
    colorcount = length(unique(TPmelt$variable))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    
    ggplot(TPmelt,aes(x=Neval,y=as.numeric(value),color=variable)) +
      geom_line() +
      geom_hline(yintercept = contact_expectation,linetype=3) +
      scale_x_continuous(breaks = protein_length*seq(0,N_TPrate,0.5),labels=seq(0,N_TPrate,0.5),expand = c(0.01,0)) +
      scale_y_continuous(expand = c(0.01,0)) +
      coord_cartesian(ylim=c(0,1)) +
      geom_vline(xintercept = protein_length*c(0.5,1),linetype=2) +
      labs(title=paste0(prefix,", lindist > ",lindist,", ",dist_cols[d]," < ",dist_cutoff,"A"),
           x="(# top contacts) / (protein length L)", y="precision (distance < 8A)", color="") +
      scale_color_manual(values=getPalette(colorcount),breaks = measured.vars, labels = label_cols0)
    ggsave(file=paste0(dataset_dir,"results/tertiary_contacts/",prefix,"TPrate_lindist",lindist,"_",dist_cols[d],"_",dist_cutoff,"A.pdf"),width=5,height =3)
    
    
    
    ##### plot eCDF
    label_cols0 = c(label_cols0,"all position pairs")
    measured.vars = c(measured.vars,"xall")
    ggplot(eCDF,aes(x=dist,color=id)) +
      stat_ecdf() +
      scale_color_manual(values=getPalette(colorcount+1),breaks = measured.vars, labels = label_cols0) +
      scale_x_continuous(breaks=seq(0,50,5),expand = c(0.01,0)) +
      scale_y_continuous(expand = c(0.01,0)) +
      geom_vline(xintercept = dist_cutoff,linetype=3) +
      labs(title=paste0(prefix,", lindist > ",lindist, ", top ",round(N_contactmap*protein_length)," contacts"),x = dist_cols[d],y = "cumulative propability",color="")
    ggsave(file=paste0(dataset_dir,"results/tertiary_contacts/",prefix,"eCDF_",dist_cols[d],"_lindist",lindist,".pdf"),width=5,height =3)
    
    
    ###### plot contactmaps
    theme_set(theme_classic())
    all_pred = measured.vars[1:(length(measured.vars)-1)]
    for (idx in seq_along(all_pred)) {
      setorderv(data2analyse,all_pred[idx],order=-1,na.last=T)
      if (modus == "cis") {
        data = data2analyse[Pos1<Pos2-lindist ][1:round(N_contactmap*protein_length),.(x=Pos1,y=Pos2,
                                                                                       shape = dist < dist_cutoff,
                                                                                       size = .SD,type=all_pred[idx]),
                                                .SDcols = all_pred[idx]]
        data = rbind(data,data2analyse[Pos1<Pos2][1:round(N_contactmap*protein_length),.(x=Pos2,y=Pos1,
                                                                                         shape = dist < dist_cutoff,
                                                                                         size = .SD,type=all_pred[idx]),
                                                  .SDcols = all_pred[idx]])
      } else if (modus == "trans") {
        data = data2analyse[1:round(N_contactmap*protein_length),.(x=Pos1,y=Pos2,
                                                                   shape = dist < dist_cutoff,
                                                                   size = .SD,type=all_pred[idx]),
                            .SDcols = all_pred[idx]]
      }
      
      data[is.na(shape),shape:=F]
      
      thisplot = ggplot()
      if (modus == "cis") {
        thisplot = thisplot + 
          geom_raster(data=data2analyse[Pos1<Pos2],aes(x=Pos1,y=Pos2,fill=dist<dist_cutoff)) +
          geom_raster(data=data2analyse[Pos1<Pos2],aes(x=Pos2,y=Pos1,fill=dist<dist_cutoff)) +
          geom_point(data=data[x<y],aes(x=x,y=y,shape=shape,size=size),color="red") +
          geom_point(data=data[x>y],aes(x=x,y=y,shape=shape,size=size),color="red") +
          geom_abline(slope=-1,intercept=0,linetype=2) +
          geom_abline(slope=-1,intercept=-lindist-.5,linetype=2,color="grey50")
      } else if (modus == "trans") {
        thisplot = thisplot + 
          geom_raster(data=data2analyse,aes(x=Pos1,y=Pos2,fill=dist<dist_cutoff)) +
          geom_point(data=data,aes(x=x,y=y,shape=shape,size=size),color="red")
      }
      if (!is.na(secondary_structure)) { #plot secondary structure above contact map
        secondary_structure[,ss_idx := rleid(ss)]
        
        if (nrow(secondary_structure[ss=="C"])>0) {
          thisplot = thisplot + 
            geom_segment(data = secondary_structure[ss=="C",.(x = min(position)-0.5,xend = max(position)+0.5,y=0,yend=0),ss_idx],
                         aes(x=x,xend=xend,y=y,yend=yend),size=0.3)
        }
        if (nrow(secondary_structure[ss=="E"])>0) {
          thisplot = thisplot + 
            geom_segment(data = ss_1pga[ss=="E",.(x = min(position)-0.5,xend = max(position)+0.5,y=0,yend=0),ss_idx],
                         aes(x=x,xend=xend,y=y,yend=yend),size=1)
        }
        if (nrow(secondary_structure[ss=="H"])>0) {
          thisplot = thisplot + 
            geom_line(data = ss_1pga[ss=="H",.(x = seq(min(position)-0.5,max(position)+0.5,0.1),
                                               y=sin(2*pi/((max(position)-min(position)+1)/round((max(position)-min(position)+1)/3.6))*seq(0,max(position)-min(position)+1,0.1))),
                                     ss_idx],
                      aes(x,y),size=0.3)
        }
      }
      thisplot = thisplot + 
        scale_fill_manual(values = c("white","grey")) +
        scale_shape_manual(breaks=c(TRUE,FALSE),values = c(1,19)) +
        scale_size_continuous(range=c(0.1,4)) +
        coord_cartesian(xlim = c(data2analyse[,min(Pos1)-0.5],data2analyse[,max(Pos1)+0.5]), ylim = c(data2analyse[,min(Pos2)-3.5],data2analyse[,max(Pos2)+0.5])) +
        scale_x_continuous(breaks=seq(5,data2analyse[,max(Pos1)],5),expand = c(0.01,0)) +
        scale_y_reverse(breaks=seq(5,data2analyse[,max(Pos2)],5),expand = c(0.01,0)) +
        labs(x="Position",y="Position",fill="distance",
             title=paste0(prefix,", ",all_pred[idx]),
             shape=paste0("<",dist_cutoff,"A"))
      
      ggsave(paste0(dataset_dir,"results/tertiary_contacts/",prefix,"contactmap_lindist",lindist,"_", dist_cols[d],"_",dist_cutoff,"A_",all_pred[idx],"_top",round(N_contactmap*protein_length),"contacts.pdf"),
             width=6,height=5)
    }
  }
}


####################################################################
##### minimal # edges < dist_cutoff between predicted contactcs ####
####################################################################

evaluate_contatcs_minimaledges = function(PWI,
                                          contactmap,
                                          dataset_dir,
                                          prefix,
                                          N_contacts = 1,
                                          lindist = 5,
                                          dist_cutoff= 8) {
  
  ### variables 
  # PWI: pairwise interaction score data.table; except for position, wt_aa and #epistatic variant columns, this should only contain the scores that should be evaluated
  # contactmap: contactmap data.table (pairwise positions) from a PDB file; can have multiple distance columns to be evaluated
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # N_contacts: number of top contacts * protein length to evaluate
  # lindist: only use position pairs with a sequence separation greater than this (only for cis libraries)
  # dist_cutoff: distances cutoff to call a direct contact (in Angstroms)
  
  #this doesn't work for trans interactions
  
  eval_cols = setdiff(names(PWI),c("Pos1","Pos2","WT_AA1","WT_AA2","NposE","NnegE"))
  dist_cols = setdiff(names(contactmap),c("Pos1","Pos2","WT_AA1","WT_AA2"))
  
  PWI2 = merge(PWI,contactmap,by=c("Pos1","Pos2"))
  
  #define protein length
  protein_length = PWI2[,diff(range(Pos1))+1] 
  
  ## calculate neighbouring matrix  using floyd-warshall algorithm
  floyd_warshall_algorithm = function(DM) { #to calculate minimal distances
    diag(DM) = 0
    for (k in 1:nrow(DM)) {
      for (i in 1:nrow(DM)) {
        for (j in 1:nrow(DM)) { 
          DM[i,j] = min(c(DM[i,j],DM[i,k]+DM[k,j]))
        } 
      } 
    }
    return(DM)
  }
  
  for (D in dist_cols) {
    # from minimal distances
    dist_dt = PWI2[,.SD,,.SDcols = c("Pos1","Pos2",D)]
    setkey(dist_dt,Pos1,Pos2)
    dist_dt = dist_dt[.(rep(unique(dist_dt$Pos1),length(unique(dist_dt$Pos1))),
                        rep(unique(dist_dt$Pos1),each=length(unique(dist_dt$Pos1)))),
                      get(D),.EACHI]
    names(dist_dt)[3] = D
    NM = matrix(dist_dt[,get(D)],nrow = length(unique(dist_dt$Pos1)), ncol = length(unique(dist_dt$Pos1)))
    NM[NM < dist_cutoff] = 1
    NM[NM >= dist_cutoff] = Inf
    dist_dt$min_edges = c(floyd_warshall_algorithm(NM))
    
    ###  edge map
    ggplot(dist_dt,aes(x=Pos1,Pos2,fill=factor(min_edges))) + 
      geom_raster() +
      scale_fill_brewer(palette = "Set1") +
      labs(fill="edges",title=D) +
      theme_classic()
    ggsave(paste0(dataset_dir,"results/tertiary_contacts/",prefix,"edgemap_",D,"_",dist_cutoff,"A.pdf"),width=5,height=4)
    
    ### barplots for min # of edges
    edge_DT = merge(PWI2[Pos1<Pos2-lindist],dist_dt[,.(Pos1,Pos2,min_edges)],by=c("Pos1","Pos2"))
    
    ### number edges for top pairs enriched versus rest, no linear chain restrictions/ ALL possible pairs
    edge_bar = edge_DT[,.N,min_edges][,':=' (freq=N/sum(N),type="all")]
    Npairs = round(N_contacts*protein_length)
    for (E in eval_cols) {
      setorderv(edge_DT,E,order = -1,na.last = T)
      edge_bar = rbind(edge_bar,edge_DT[get(E) >= edge_DT[Npairs,get(E)],.N,
                                        by=.(min_edges)][,':=' (freq=N/sum(N),type=E)])
    }
    
    ggplot(edge_bar,aes(x=type,y=freq,fill=as.ordered(min_edges))) +
      geom_col(position = position_stack(reverse = TRUE),width=0.75) +
      scale_fill_brewer(direction=-1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0("top ",Npairs," pairs > ",lindist,"aa"),x="",y="fraction",fill=paste0("# edges < ",dist_cutoff,"A"))
    ggsave(file=paste0(dataset_dir,"results/tertiary_contacts/",prefix,"minNedges_top",Npairs,"_",D,"_",dist_cutoff,"A_lindist",lindist,".pdf"),width=5,height=5)
  }
}



####################################################
##### plot interaction scores versus distances #####
####################################################
score_vs_distance_scatter = function(contacts,
                                        contactmap,
                                        dataset_dir,
                                        prefix,
                                        modus = "cis",
                                        lindist = 5,
                                        dist_cutoff= 8) {
  
  ### variables 
  # contacts: pairwise interaction score data.table; except for position, wt_aa and #epistatic variant columns, this should only contain the scores that should be evaluated
  # contactmap: contactmap data.table (pairwise positions) from a PDB file; can have multiple distance columns to be evaluated
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # lindist: only use position pairs with a sequence separation greater than this (only for cis libraries)
  # dist_cutoff: distances cutoff to call a direct contact (in Angstroms)
  # modus: "cis" (single protein) or "trans" (protein-protein interaction)
  
  
  require(data.table)
  require(gtable)
  require(cowplot)
  # require(mixtools)
  theme_set(theme_classic(base_size = 9))
  
  if (modus == "cis") {
    contacts = contacts[Pos1<Pos2-lindist]
    contactmap = contactmap[Pos1<Pos2-lindist]
  }
  
  eval_cols = setdiff(names(contacts),c("Pos1","Pos2","WT_AA1","WT_AA2","NposE","NnegE"))
  dist_cols = setdiff(names(contactmap),c("Pos1","Pos2","WT_AA1","WT_AA2"))
  
  for (i in seq_along(eval_cols)) {
    data = contacts[,.(Pos1,Pos2,.SD),,.SDcols = eval_cols[i]]
    names(data)[3] = "predictor"
    
    for (d  in seq_along(dist_cols)) {
      
      data2 = merge(data,contactmap[,.(Pos1,Pos2,.SD),,.SDcols = dist_cols[d]])  
      names(data2)[4] = "distance"
      
      data2[,bin_dist := findInterval(distance,seq(0,max(distance),dist_cutoff))]
      
      # mixEM = mixtools::normalmixEM(data2[!is.na(predictor),predictor],k=2,mean.constr = c(NA,">a"))
      # xrange = data2[!is.na(predictor),seq(min(predictor),max(predictor),length.out = 100)]
      # p1 = ggplot() + 
      #   geom_density(data = data,aes(predictor)) +
      #   geom_line(data = data.frame(x = xrange,
      #                               d = dnorm(x = xrange,mean = mixEM$mu[which.min(mixEM$mu)],sd = mixEM$sigma[which.min(mixEM$mu)])),
      #             aes(x,y=d),color="grey") +
      #   geom_line(data = data.frame(x = xrange,
      #                               d = dnorm(x = xrange,mean = mixEM$mu[which.max(mixEM$mu)],sd = mixEM$sigma[which.max(mixEM$mu)])),
      #             aes(x,y=d),color="red") +
      #   scale_x_continuous(expand = c(0.01,0)) +
      #   labs(x = eval_cols[i])
      # 
      # ref = 0.9
      # ggplot(data2[!is.na(predictor)]) + 
      #   geom_point(aes(y=predictor,x=distance,color = mixEM$posterior[,which.max(mixEM$mu)])) +
      #   geom_boxplot(aes(y = predictor,x=distance,group = bin_dist),alpha=0,outlier.alpha = 0) +
      #   geom_vline(xintercept = dist_cutoff,linetype=2)+
      #   coord_flip() +
      #   scale_color_gradient(low="grey",high="red") +
      #   # scale_color_gradientn(colours = c("grey","grey","red","red","red"),
      #   # values=c(0,ref-0.05,ref,ref+0.05,1)) +
      #   scale_y_continuous(expand = c(0.01,0)) +
      #
      #   labs(y = eval_cols[i],x = dist_cols[d],color = "P(contact)",
      #        title = data2[!is.na(predictor)][mixEM$posterior[,which.max(mixEM$mu)]>0.9,.(.N,TP = sum(distance < dist_cutoff))][
      #          ,paste0("P > 0.9: ",TP,"/",N," TP (",round(TP/N*100),"%)")])
      # # plot_grid(plotlist = list(p1,p2),nrow=2,rel_heights = c(1,2),axis = "l")
      
      ggplot(data2[!is.na(predictor)]) +
        geom_point(aes(y=predictor,x=distance),color="grey75") +
        geom_boxplot(aes(y = predictor,x=distance,group = bin_dist),alpha=0,outlier.alpha = 0) +
        geom_vline(xintercept = dist_cutoff,linetype=2)+
        coord_flip() +
        scale_y_continuous(expand = c(0.02,0)) +
        scale_x_continuous(expand = c(0.02,0)) +
        labs(y = eval_cols[i],x = dist_cols[d],
             title = paste0("R = ",data2[,round(cor(predictor,distance),digits=2)]))
      ggsave(file = paste0(dataset_dir,"results/tertiary_contacts/",prefix,"scatter_",eval_cols[i],"_",dist_cols[d],".pdf"),width=4,height=4)
      
    }
  }
}
