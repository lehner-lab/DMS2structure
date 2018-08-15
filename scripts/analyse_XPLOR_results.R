#analyse results from XPLOR simulations

analyse_XPLOR_results = function(XPLOR_dir,
                                 filepattern = ".RData",
                                 contactmap,
                                 best_models = 0.05,
                                 eval_L = 1,
                                 draw_contactmaps = FALSE) {
  
  ### variables
  # XPLOR_dir: directory where the data to be analysed is stored, will create a subfolder "results" for all plots and tables
  # filepattern: in XPLOR_dir, analyse all files that match the filepattern
  # contactmap: pair-wise distance data.table from reference structure to compare to (draw contactmaps, precision changes over stages)
  # best_models: fraction of top energy models to evaluate in terms of RMSD and TMscore
  # eval_L: for some plots, for which simulations (eval_L * protein length contacts used for distance restraints) to show top models for
  # draw_contactmaps: if TRUE, will draw a contactmap with initially predicted restraints and reweighted restraints at stage 3 for every file
  
  require(data.table)
  require(ggplot2)
  require(cowplot)
  theme_set(theme_minimal())
  
  filenames = dir(path = XPLOR_dir,pattern = filepattern)
  
  results=list()
  for (i in seq_along(filenames)) {
    load(file=paste0(XPLOR_dir,filenames[i]))
    results[[i]] = varlist
    results[[i]]$simulation = paste0(results[[i]]$predictor,"_L",results[[i]]$L)
    if (i ==1 ) {
      ID_table = data.table(idx = i,predictor = results[[i]]$predictor, L = results[[i]]$L)
    } else {
      ID_table = rbind(ID_table,data.table(idx = i,predictor = results[[i]]$predictor, L = results[[i]]$L))  
    }
  }
  
  system(command = paste0("mkdir ",XPLOR_dir,"results/"))
  
  ############ plot results ############
  
  #############################################################
  ##### RMSD as function of structural model total energy #####
  #############################################################
  DT_RMSD = data.table(predictor = "",L = as.numeric(NA),total = as.numeric(NA),RMSD = as.numeric(NA),stage = as.numeric(NA),cluster = as.numeric(NA),shape = NA)
  for (i in seq_along(filenames)) {
    for (stage in 1:3) {
      DT_RMSD = rbind(DT_RMSD,results[[i]]$energy_XPLOR[[stage]][cluster!="avg",.(predictor = results[[i]]$predictor,L = results[[i]]$L,
                                                                                  total,RMSD = TM_RMSD,stage,cluster, 
                                                                                  shape = total < results[[i]]$energy_XPLOR[[stage]][min_cluster_top10==T,max(total)])])  
    }
  }
  DT_RMSD = DT_RMSD[!is.na(predictor)]
  
  unique_predictors = unique(ID_table$predictor)
  for (P in unique_predictors) {
    ggplot(DT_RMSD[predictor == P],aes(total,RMSD,color=cluster, shape = shape,alpha = shape)) +
      geom_point() +
      scale_shape_manual(values = c(1,16)) +
      scale_alpha_manual(values = c(.25,1)) +
      coord_cartesian(xlim = quantile(DT_RMSD[predictor == P,total],c(0,.95)),ylim = c(0,DT_RMSD[predictor == P,max(RMSD)])) +
      labs(x="total XPLOR energy",y="accuracy",color="cluster",shape = "top10percent") +
      facet_grid(L ~ stage) 
    ggsave(paste0(XPLOR_dir,"results/RMSD_vs_totalenergy_",P,".pdf"),width=8,height=5)
  }
  
  
  #######################################################
  ##### precision of predicted contacts over stages #####
  #######################################################
  
  DT_precision = data.table(predictor = "",L = as.numeric(NA),topN = as.numeric(NA),stage = as.numeric(NA),precision = as.numeric(NA))
  for (i in seq_along(filenames)) {
    if (results[[i]]$predictor != "control") {
      
      # helper = merge(results[[i]]$NOE_DT[[1]][type == "contact",-grep("restraint_NR",names(results[[i]]$NOE_DT[[1]])),with=F],contactmap,by=c("Pos1","Pos2"))
      helper = merge(results[[i]]$NOE_DT[[1]][type == "contact",.(Pos1,Pos2,weight)],contactmap,by=c("Pos1","Pos2"))
      setorderv(helper,"weight",order=-1,na.last=T)
      DT_precision = rbind(DT_precision,data.table(predictor = results[[i]]$predictor,L = results[[i]]$L,stage = "initial weighting",
                                                   topN = 1:nrow(helper),precision = helper[,cumsum(scHAmin<8)/(1:.N)]))
      
      # helper = merge(results[[i]]$NOE_DT[[3]][type == "contact" & weight > 0.1,-grep("restraint_NR",names(results[[i]]$NOE_DT[[3]])),with=F],contactmap,by=c("Pos1","Pos2"))
      helper = merge(results[[i]]$NOE_DT[[3]][type == "contact" & weight0 > 0.1,.(Pos1,Pos2,weight0)],contactmap,by=c("Pos1","Pos2"))
      setorderv(helper,"weight0",order=-1,na.last=T)
      DT_precision = rbind(DT_precision,data.table(predictor = results[[i]]$predictor,L = results[[i]]$L,stage = "final weighting",
                                                   topN = 1:nrow(helper),precision = helper[,cumsum(scHAmin<8)/(1:.N)]))
    }
  }
  
  DT_precision = DT_precision[!is.na(predictor)]
  
  unique_predictors = unique(ID_table$predictor)
  for (P in unique_predictors) {
    if (P != "control") {
      ggplot(DT_precision[predictor == P],aes(topN/results[[i]]$protein_length,precision,color=factor(L), linetype = stage)) +
        geom_line() +
        scale_color_brewer(palette = "Set1") +
        coord_cartesian(ylim = c(0,1)) +
        labs(x="top predicted contacts",y="precision",color="# contacts used",linetype = "stage") +
        facet_wrap(~L)
      ggsave(paste0(XPLOR_dir,"results/precision_",P,".pdf"),width=8,height=5)  
    }
  }
  
  #############################################################
  ##### contact map with predictions stage1 versus stage3 #####
  #############################################################
  if (draw_contactmaps) {
    for (i in seq_along(filenames)) {
      if (!grepl("control",filenames[i])) {
        #contact map with predictions stage1 versus stage3
        max_HA = contactmap[,max(scHAmin,na.rm=T)]
        
        NOE_DT2 = merge(results[[i]]$NOE_DT[[2]][weight0>0.1 &type == "contact",.(Pos1,Pos2,weight0)],
                        contactmap,by=c("Pos1","Pos2"))
        NOE_DT3 = merge(results[[i]]$NOE_DT[[3]][weight>0.1 & type == "contact",.(Pos1,Pos2,weight)],
                        contactmap,by=c("Pos1","Pos2"))
        
        PLOT = ggplot() + 
          geom_raster(data=contactmap[Pos1<Pos2],aes(x=Pos1,y=Pos2,fill=scHAmin<8)) +
          geom_raster(data=contactmap[Pos1<Pos2],aes(x=Pos2+1,y=Pos1-1,fill=scHAmin<8)) +
          scale_fill_manual(breaks = c(FALSE,TRUE),values=c("white", "grey75")) +
          ## predicted contacts
          geom_point(data=NOE_DT2,aes(x=Pos1,y=Pos2,shape=scHAmin<8,size=weight0),color="red") +
          geom_point(data=NOE_DT3,aes(x=Pos2+1,y=Pos1-1,shape=scHAmin<8,size=weight),color="orange") +
          ## predicted secondary structure
          geom_point(data=results[[i]]$DIHE_DT[[2]][weight0>0.33,w := mean(weight0),.(position,ss)],aes(x=position,y=position,size=w,shape=ss),color="black") +
          geom_point(data=results[[i]]$DIHE_DT[[3]][weight>0.33,w := mean(weight),.(position,ss)],aes(x=position+1,y=position-1,size=w,shape=ss),color="black") +
          
          scale_shape_manual(breaks=c(FALSE,TRUE,"alpha","beta"),values = c(1,19,15,17)) +
          scale_size_continuous(range=c(0.1,2.5)) +
          # scale_color_brewer(palette = "Set1") +
          coord_cartesian(xlim = c(0.5,results[[i]]$protein_length+1.5), ylim = c(0.5,results[[i]]$protein_length+0.5)) + 
          geom_abline(slope=-1,intercept=1,linetype=2) +
          geom_abline(slope=-1,intercept=-5.5,linetype=2,color="grey50") +
          geom_abline(slope=-1,intercept=7.5,linetype=2,color="grey50") +
          scale_x_continuous(breaks=seq(5,results[[i]]$protein_length,5),expand = c(0.005,0)) +
          scale_y_reverse(breaks=seq(5,results[[i]]$protein_length,5),expand = c(0.005,0)) +
          labs(x="position",y="position",fill="distance",
               title=paste0(results[[i]]$simulation,", initial(red) vs. final(orange)"),
               shape="< 8A")
        
        #predicted beta sheet hbonds
        if (results[[i]]$NOE_DT[[1]][type == "beta sheet hbond",.N > 0]) {
          PLOT = PLOT + geom_point(data=results[[i]]$NOE_DT[[1]][type == "beta sheet hbond" & Pos1 < Pos2],aes(x=Pos1,y=Pos2,size=1/2),color="black",alpha=0.5,shape=18) +
            geom_point(data=results[[i]]$NOE_DT[[1]][type == "beta sheet hbond" & Pos2 < Pos1],aes(x=Pos2,y=Pos1,size=1/2),color="black",alpha=0.5,shape=18) +
            geom_point(data=results[[i]]$NOE_DT[[3]][weight>0.1 & type == "beta sheet hbond" & Pos1 < Pos2],aes(x=Pos2+1,y=Pos1-1,size=weight/2),color="black",alpha=0.5,shape=18) +
            geom_point(data=results[[i]]$NOE_DT[[3]][weight>0.1 & type == "beta sheet hbond" & Pos2 < Pos1],aes(x=Pos1+1,y=Pos2-1,size=weight/2),color="black",alpha=0.5,shape=18)
          
          
          if (results[[i]]$predictor != "control") {
            PLOT = PLOT + geom_point(data=results[[i]]$NOE_DT[[1]][type == "beta sheet hbond" & Pos1_opt2 < Pos2_opt2],aes(x=Pos1_opt2,y=Pos2_opt2,size=1/2),color="black",alpha=0.5,shape=18) +
              geom_point(data=results[[i]]$NOE_DT[[1]][type == "beta sheet hbond" & Pos2_opt2 < Pos1_opt2],aes(x=Pos2_opt2,y=Pos1_opt2,size=1/2),color="black",alpha=0.5,shape=18) +
              geom_point(data=results[[i]]$NOE_DT[[3]][weight>0.1 & type == "beta sheet hbond" & Pos1_opt2 < Pos2_opt2],aes(x=Pos2_opt2+1,y=Pos1_opt2-1,size=weight/2),color="black",alpha=0.5,shape=18) +
              geom_point(data=results[[i]]$NOE_DT[[3]][weight>0.1 & type == "beta sheet hbond" & Pos2_opt2 < Pos1_opt2],aes(x=Pos1_opt2+1,y=Pos2_opt2-1,size=weight/2),color="black",alpha=0.5,shape=18)
          }
        }
        PLOT
        ggsave(paste0(XPLOR_dir,"results/contactmap_reweighting_",results[[i]]$simulation,".pdf"),width=6,height=5)
      }
    }
  }
  
  
  
  #########################################
  ##### accuracy of structural models #####
  #########################################
  
  DT_structmodels = data.table(predictor = "",DC = NA,L = as.numeric(NA),total= as.numeric(NA),pdb_id = "",
                               TM_score = as.numeric(NA),RMSD = as.numeric(NA),stage = as.numeric(NA))
  for (i in seq_along(filenames)) {
    for (s in 1:3) {
      DT_structmodels = rbind(DT_structmodels,
                              results[[i]]$energy_XPLOR[[s]][order(total)][1:round(.N*best_models),.(predictor = gsub("_DC","",results[[i]]$predictor),
                                                                                                     DC = grepl("_DC",results[[i]]$predictor),
                                                                                                     L = results[[i]]$L,
                                                                                                     total,pdb_id,TM_score,RMSD=TM_RMSD,stage=s)])
    }
  }
  DT_structmodels = DT_structmodels[!is.na(DC)]
  DT_structmodels[stage==3,RMSD_mean := mean(RMSD),.(predictor,L,DC)]
  DT_structmodels[stage==3,TM_score_mean := mean(TM_score),.(predictor,L,DC)]
  
  setkey(DT_structmodels,predictor,L,DC,total)
  write.table(x = DT_structmodels[stage==3,.(predictor,DC,L,total,pdb_id,TM_score,TM_score_mean,RMSD,RMSD_mean)],
              file = paste0(XPLOR_dir,"results/top_models.txt"),sep = "\t",quote = FALSE,row.names = FALSE)
  
  ################
  # models as function of predictor, deep contact and number of restraints used

  DT_structmodels[,xlabel := paste0(DC,"_",L,collapse=""),.(DC,L)]
  
  ### TM score
  ggplot(DT_structmodels[stage==3],aes(x=xlabel,y=TM_score,fill=L,color = DC)) +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_distiller(palette = "Blues",direction = 1) +
    scale_color_manual(values = c("orange","red")) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    coord_cartesian(ylim = c(0,1)) +
    facet_wrap( ~ predictor) +
    labs(y="TM score",
         title=paste0(results[[1]]$protein,", ",best_models*100,"% best models"),
         color = "deepcontact",fill="#restraints/L")
  ggsave(paste0(XPLOR_dir,"results/TMscore_top",best_models*100,"per_models.pdf"),width=13,height=13)
  
  ### RMSD
  ggplot(DT_structmodels[stage==3],aes(x=xlabel,y=RMSD,fill=L,color = DC)) +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_distiller(palette = "Blues",direction = 1) +
    scale_color_manual(values = c("orange","red")) +
    scale_y_log10(breaks = c(1,2,3,4,5,7.5,10)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    facet_wrap( ~ predictor) +
    labs(y="RMSD",
         title=paste0(results[[1]]$protein,", ",best_models*100,"% best models"),
         color = "deepcontact",fill="#restraints/L")
  ggsave(paste0(XPLOR_dir,"results/RMSD_top",best_models*100,"per_models.pdf"),width=13,height=13)
  
  
  
  ################
  # models as function of predictor and deep contact but aggregated across number of restraints used
  
  require(RColorBrewer)
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  DT_structmodels[,xlabel := NULL]
  DT_structmodels[,xlabel := paste0(predictor,"_",DC,collapse=""),.(predictor,DC)]
  DT_structmodels[,xlabel := factor(xlabel,levels = unique(DT_structmodels[,.(xlabel,predictor)])[order(predictor)]$xlabel)]
  
  ### TMscore
  ggplot(DT_structmodels[stage==3 & predictor != "control"],aes(x=xlabel,y=TM_score,color=predictor,fill = DC)) +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_manual(values = c("white","grey75")) +
    scale_color_manual(values=getPalette(DT_structmodels[stage==3 & predictor != "control",length(unique(predictor))])) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(y="TM score",
         title=paste0(results[[1]]$protein,", ",best_models*100,"% best models, all L"),
         fill = "deepcontact",color="predictor")
  ggsave(paste0(XPLOR_dir,"results/TMscore_top",best_models*100,"per_models_aggregated.pdf"),width=8,height=6)
  
  
  ### RMSD
  ggplot(DT_structmodels[stage==3 & predictor != "control"],aes(x=xlabel,y=RMSD,color=predictor,fill = DC)) +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_manual(values = c("white","grey75")) +
    scale_color_manual(values=getPalette(DT_structmodels[stage==3 & predictor != "control",length(unique(predictor))])) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    coord_cartesian(ylim = c(0,DT_structmodels[stage==3 & predictor != "control",max(RMSD)])) +
    labs(y="RMSD",
         title=paste0(results[[1]]$protein,", ",best_models*100,"% best models, all L"),
         fill = "deepcontact",color="predictor")
  ggsave(paste0(XPLOR_dir,"results/RMSD_top",best_models*100,"per_models_aggregated.pdf"),width=8,height=6)
  
  

  ################
  # models as function of predictor and deep contact at L = eval_L
  
  require(RColorBrewer)
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  DT_structmodels[,xlabel := NULL]
  DT_structmodels[,xlabel := paste0(predictor,"_",DC,collapse=""),.(predictor,DC)]
  DT_structmodels[,xlabel := factor(xlabel,levels = unique(DT_structmodels[,.(xlabel,predictor)])[order(predictor)]$xlabel)]
  
  ### TMscore
  ggplot(DT_structmodels[stage==3 & L == eval_L],aes(x=xlabel,y=TM_score,color=predictor,fill = DC)) +
    geom_boxplot(outlier.shape=NA) + 
    scale_fill_manual(values = c("white","grey75")) +
    scale_color_manual(values=getPalette(DT_structmodels[stage==3 & L == eval_L,length(unique(predictor))])) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(y="TM score",
         title=paste0(results[[1]]$protein,", ",best_models*100,"% best models, L = ",eval_L),
         fill = "deepcontact",color="predictor")
  ggsave(paste0(XPLOR_dir,"results/TMscore_top",best_models*100,"per_models_L",eval_L,".pdf"),width=8,height=6)
  
  
  ### RMSD
  ggplot(DT_structmodels[stage==3 & L == eval_L],aes(x=xlabel,y=RMSD,color=predictor,fill = DC)) +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_manual(values = c("white","grey75")) +
    scale_color_manual(values=getPalette(DT_structmodels[stage==3 & L == eval_L,length(unique(predictor))])) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    coord_cartesian(ylim = c(0,DT_structmodels[stage==3 & L == eval_L,max(RMSD)])) +
    labs(y="RMSD",
         title=paste0(results[[1]]$protein,", ",best_models*100,"% best models, L = ",eval_L),
         fill = "deepcontact",color="predictor")
  ggsave(paste0(XPLOR_dir,"results/RMSD_top",best_models*100,"per_models_L",eval_L,".pdf"),width=8,height=6)
  
  
  ################
  # evolution of model accuracy across modeling stages
  
  DT_structmodels[,predictor_DC := paste0(predictor,ifelse(DC,"_DC",""),collapse=""),.(predictor,DC)]
  DT_structmodels[,xlabel := NULL]
  DT_structmodels[,xlabel := paste0("L",L,"_stage",stage,collapse=""),.(L,stage)]
  DT_structmodels[,xlabel := factor(xlabel,levels = unique(DT_structmodels[,.(xlabel,L,stage)])[order(L,stage)]$xlabel)]
  
  #### TMscore
  ggplot(DT_structmodels[predictor != "control"],aes(x=xlabel,y=TM_score,fill=factor(stage),color=factor(L))) +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_brewer(palette="Blues") +
    scale_color_manual(values=getPalette(DT_structmodels[,length(unique(L))])) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    facet_wrap(~ predictor_DC) +
    labs(fill = "stage",color = "#restraints")
  ggsave(paste0(XPLOR_dir,"results/TMscore_stage_improvement_top",best_models*100,"per_models_.pdf"),width=12,height=12)
  
  #### RMSD
  ggplot(DT_structmodels[predictor != "control"],aes(x=xlabel,y=RMSD,fill=factor(stage),color=factor(L))) +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_brewer(palette="Blues") +
    scale_color_manual(values=getPalette(DT_structmodels[,length(unique(L))])) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    coord_cartesian(ylim = c(0,quantile(DT_structmodels[predictor != "control",RMSD],.99))) +
    facet_wrap(~ predictor_DC) +
    labs(fill = "stage",color = "#restraints")
  ggsave(paste0(XPLOR_dir,"results/RMSD_stage_improvement_top",best_models*100,"per_models_.pdf"),width=12,height=12)
}
