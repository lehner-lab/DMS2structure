######################################################################
### use DeepContact convolutional neural network (Liu et al. 2107) ###
### to transform interaction score matrices; this basic architecture #
### only uses the 2d interaction score matrix as input and outputs ###
### a transformed interaction score matrix ###########################
######################################################################


deepcontact_transform_basic2d = function(PWI,
                                 dataset_dir,
                                 output_filename = "DMS_PWI_deepcontact.txt",
                                 prefix = "",
                                 modus = "cis",
                                 plot_results = T,
                                 deepcontact_dir,
                                 normalize = T) {
  
  ### variables 
  # PWI: pairwise interaction score data.table
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # output_filename: filename to write datatable to in dataset_dir/processed_data/
  # prefix:to be added to results files (in case of running diff. versions of data from same dataset etc)
  # modus: "cis" (single protein) or "trans" (protein-protein interaction)
  # plot results: if TRUE, will output untransformed versus transformed data into dataset_dir/results/deepcontact/
  # normalize: normalize interaction scores to range from 0 to 1 (as do CCMpred scores, which are normally used with DeepContact)
  # deepcontact_dir: directory where the DeepContact repository was cloned to (use "git clone https://github.com/largelymfs/deepcontact")
  #####   NOTE that you'll also have to install the python dependencies required for DeepContact (see github page point 4; points 1-3 can be ignored for this basic "CCMpred only" architecture) 
  
  
  require(gdata)
  require(cowplot)
  
  #setup deep contact
  Sys.setenv(PYTHONPATH = deepcontact_dir)
  if (length(grep("anaconda",Sys.getenv("PATH")))==0) {
    Sys.setenv(PATH = paste0(Sys.getenv("PATH"),":/anaconda/envs/deepcontact-env/bin"))
  }
  
  #write input scores to files
  setkey(PWI,Pos1,Pos2)
  unique_pos = unique(c(PWI$Pos1,PWI$Pos2))
  protein_length = max(unique_pos) - min(unique_pos) + 1
  
  #set up data.table for output
  PWI_transformed = data.table(Pos1 = rep(unique_pos,protein_length),
                               Pos2 = rep(unique_pos,each=protein_length))
  
  #find interaction score columns for evaluation
  eval_cols = setdiff(names(PWI),c("Pos1","Pos2","WT_AA1","WT_AA2","NposE","NnegE"))
  
  for (i in seq_along(eval_cols)) {
    setkey(PWI,Pos1,Pos2)
    setkey(PWI_transformed,Pos1,Pos2)
    
    if (modus == "cis") {
      if (normalize == T) {#normalize
        PWI[Pos1 != Pos2,(eval_cols[i]) := (.SD - min(.SD,na.rm=T)) / (max(.SD,na.rm=T) - min(.SD,na.rm=T)),.SDcols = eval_cols[i]]
      }
      #set diagonal elements to 0
      PWI[Pos1 == Pos2,(eval_cols[i]) := 0]
    } else { #if trans, don't treat diagonal separately
      if (normalize == T) {#normalize
        PWI[,(eval_cols[i]) := (.SD - min(.SD,na.rm=T)) / (max(.SD,na.rm=T) - min(.SD,na.rm=T)),.SDcols = eval_cols[i]]
      }
    }
    #make matrix
    helper = matrix(PWI[.(rep(unique_pos,protein_length),rep(unique_pos,each=protein_length)),unlist(.SD),.SDcols=eval_cols[i]],nrow=protein_length,ncol=protein_length)

    if (modus == "cis") { #make symmetric
      lowerTriangle(helper) = upperTriangle(helper,byrow=T) 
    }
    
    #write to file
    system(command = paste0("mkdir ",deepcontact_dir,"ccmpred_only/"),wait=T)
    write.table(helper,file = paste0(deepcontact_dir,"ccmpred_only/tmp.ccmpred"),
                quote = F,sep = "\t",row.names = F,col.names = F)
    
    #transform with deep contact
    org_wd = getwd()
    setwd(deepcontact_dir)
    command_input = paste0("source activate deepcontact-env","\n","python ",deepcontact_dir,"scripts/predict_using_ccmpred.py --input_filename ",deepcontact_dir,"ccmpred_only/tmp.ccmpred --output_filename ",deepcontact_dir,"ccmpred_only/tmp.output")
    system(command = command_input,wait = T)
    setwd(org_wd)
    
    #read transformed files
    if (modus == "cis") {
      DC = as.matrix(fread(paste0(deepcontact_dir,"ccmpred_only/tmp.output")))
      DC_norm = (DC + t(DC))/2
      PWI_transformed[,paste0(eval_cols[i]) := c(DC_norm)]
    } else {
      PWI_transformed[,paste0(eval_cols[i]) := c(as.matrix(fread(paste0(deepcontact_dir,"ccmpred_only/tmp.output"))))]  
    }
    
    ### plot results
    if (plot_results) {
      p1 = ggplot() +
        geom_raster(data=PWI[Pos1<Pos2],aes(Pos1,Pos2,fill=get(eval_cols[i]))) +
        geom_raster(data=PWI[Pos1<Pos2],aes(Pos2,Pos1,fill=get(eval_cols[i]))) +
        scale_fill_gradient(low = "white",high="dodgerblue1") +
        geom_abline(slope=-1,linetype=2) +
        scale_x_continuous(breaks = seq(0,protein_length,5),expand = c(0,0)) +
        scale_y_reverse(breaks = seq(0,protein_length,5),expand = c(0,0)) +
        labs(fill = "",x = "position", y = "position",
             title = eval_cols[i])
      p2 = ggplot() +
        geom_raster(data=PWI_transformed,aes(Pos1,Pos2,fill=get(eval_cols[i]))) +
        scale_fill_gradient(low = "white",high="dodgerblue1") +
        geom_abline(slope=-1,linetype=2) +
        scale_x_continuous(breaks = seq(0,protein_length,5),expand = c(0,0)) +
        scale_y_reverse(breaks = seq(0,protein_length,5),expand = c(0,0)) +
        labs(fill = "",x = "position", y = "position",
             title = "DeepContact")
      plot_grid(plotlist = list(p1,p2),nrow=1)
      ggsave(paste0(dataset_dir,"results/deepcontact/",prefix,eval_cols[i],".pdf"),width=11,height=5)
      
    }
  }
  
  # write transformed data.table to file
  write.table(x = PWI_transformed, file = paste0(dataset_dir,"processed_data/",output_filename),
              quote = F,row.names = F, col.names = T)
  
  # return transformed data.table
  return(PWI_transformed)
}