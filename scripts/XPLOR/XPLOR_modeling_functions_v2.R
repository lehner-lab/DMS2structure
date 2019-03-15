
XPLOR_modeling = function(varlist) {
  
  require(data.table)
  
  #create working directories
  system(command = paste0("mkdir ",varlist$folder,"XPLOR_output"), wait = T)
  system(command = paste0("mkdir ",varlist$folder,"R_output"), wait = T)
  system(command = paste0("mkdir ",varlist$folder,"top_structures"), wait = T)
  
  #start xplor modeling
  for (stage in 1:3) {
    
    #write NOE for contacts and SS
    varlist = write_contact_NOE_restraints(varlist,stage)
    #write dihe restraints
    write_SS_dihe_restraints(varlist,stage)
    #write XPLOR python script
    modify_XPLOR_phython(varlist,stage)
    
    ############ run XPLOR ############
    #server-shell: start simulations
    system(paste0("/software/bl/el7.2/xplor-nih-2.47/bin/xplor -smp ",varlist$cores," -py -o ",varlist$folder,varlist$filename[stage],".out ",varlist$folder,varlist$filename[stage],".py"),
           wait = T)
    
    #check whether all structures were calculated, if not repeat
    NOS_exp = varlist$numberOfStructures[stage]
    NOS_obs = length(grep(paste0(varlist$filename[stage],"_[0-9]+.*pdb[_]*[1-9]*$"),list.files(path = paste0(varlist$folder,"XPLOR_output/"))))
    tries = 1
    while (NOS_obs < NOS_exp) {
      if (tries < 5) {
        varlist$numberOfStructures[stage] = NOS_exp - NOS_obs
        print(paste0("XPLOR fail: missing ",NOS_exp - NOS_obs," structures in stage ",stage))
        modify_XPLOR_phython(varlist,stage)
        system(paste0("/software/bl/el7.2/xplor-nih-2.47/bin/xplor -smp ",varlist$cores," -py -o ",varlist$folder,varlist$filename[stage],".out ",varlist$folder,varlist$filename[stage],".py"),
               wait = T)
        tries = tries + 1
      } else { #if this doesn't work after five tries, stop the simulations
        "something went wrong with repeated attempts to generate enough structures"
        error
      }
    }
    varlist$numberOfStructures[stage] = NOS_exp
    
    
    ############ evaluate models ############ 
    #read energies files & restraint violations, evaluate top models, reweight restraints
    varlist = read_PDBfiles_violations(varlist,stage)
  }
  
  ### save collected results
  save(file = paste0(varlist$folder,"variables_results.RData"),varlist)
  
  #save results also to one-up directory to minimize copying from server
  a=gregexpr(varlist$folder,pattern = "/")[[1]]
  save(file = paste0(paste0(strsplit(varlist$folder,split="")[[1]][1:a[length(a)-2]],collapse=""),
                     varlist$protein,varlist$predictor,"_L",varlist$L,"_variables_results.RData"),varlist)
  
  
  ## copy top 5 structure from each stage to new dir
  top5 = c(varlist$energy_XPLOR[[1]][order(total)][1:5,pdb_id],
           varlist$energy_XPLOR[[2]][order(total)][1:5,pdb_id],
           varlist$energy_XPLOR[[3]][order(total)][1:5,pdb_id])
  for (i in 1:length(top5)) {
    system(command = paste0("cp ",paste0(varlist$folder,"XPLOR_output/",top5[i])," ",varlist$folder,"top_structures/"),wait=T)  
  }
  
  
  #remove large output files from XPLOR modeling after simulations are done
  system(paste0("rm ",varlist$folder,"*.out*"))
  # and remove all structural models in the XPLOR output folder
  # system(command = paste0("rm -rf ",varlist$folder,"XPLOR_output/"))
} 

write_contact_NOE_restraints = function(varlist,stage){
  
  #this file writes a NOE restraint file for distance restraints derived from tertiary contacts, beta strands and beta sheet hbonding
  
  NOE_string = "set echo=false end"
  NOE_string[2] = "set wrnlev=0 end"
  idx=3
  restraint_idx = 1
  varlist$NOE_DT[[stage]] = cbind(varlist$NOE_DT[[stage]],data.table(restraint_NR = as.numeric(NA)))
  #### predicted contacts
  if (varlist$NOE_DT[[stage]][type=="contact" & weight > 0.1,.N]>0) {
    NOE_string[idx] = "!contacts cb pairing"
    idx=idx+1
    for (i in varlist$NOE_DT[[stage]][,which(type=="contact"& weight > 0.1)]) {
      NOE_string[idx] = sprintf("WEIGht %.3f",varlist$NOE_DT[[stage]][i,weight])
      NOE_string[idx+1] = sprintf("ASSIgn (resid %i and name %s)(resid %i and name %s) %.1f %.1f %.1f",
                                  varlist$NOE_DT[[stage]][i,Pos1],varlist$NOE_DT[[stage]][i,atom1],
                                  varlist$NOE_DT[[stage]][i,Pos2],varlist$NOE_DT[[stage]][i,atom2],
                                  varlist$NOE_DT[[stage]][i,dist],varlist$NOE_DT[[stage]][i,lower_dist],varlist$NOE_DT[[stage]][i,upper_dist])
      idx=idx+2
      varlist$NOE_DT[[stage]][i,restraint_NR := restraint_idx]
      restraint_idx = restraint_idx+1
    }
  }
  
  #### intra betastrand O-O contacts
  if (varlist$NOE_DT[[stage]][type=="beta strand" & weight > 0.1,.N]>0) {
    NOE_string[idx] = "!beta strand O-O pairing"
    idx=idx+1
    for (i in varlist$NOE_DT[[stage]][,which(type=="beta strand" & weight > 0.1)]) {
      NOE_string[idx] = sprintf("WEIGht %.3f",varlist$NOE_DT[[stage]][i,weight])
      NOE_string[idx+1] = sprintf("ASSIgn (resid %i and name %s)(resid %i and name %s) %.1f %.1f %.1f",
                                  varlist$NOE_DT[[stage]][i,Pos1],varlist$NOE_DT[[stage]][i,atom1],
                                  varlist$NOE_DT[[stage]][i,Pos2],varlist$NOE_DT[[stage]][i,atom2],
                                  varlist$NOE_DT[[stage]][i,dist],varlist$NOE_DT[[stage]][i,lower_dist],varlist$NOE_DT[[stage]][i,upper_dist])
      idx=idx+2
      varlist$NOE_DT[[stage]][i,restraint_NR := restraint_idx]
      restraint_idx = restraint_idx+1
    }
  }
  
  #### betasheet hydrogenbonding
  if (varlist$NOE_DT[[stage]][type=="beta sheet hbond" & weight > 0.1,.N]>0) {
    NOE_string[idx] = "!beta sheet hydrogen bonds"
    idx = idx + 1
    for (i in varlist$NOE_DT[[stage]][,which(type=="beta sheet hbond" & weight > 0.1)]) {
      NOE_string[idx] = sprintf("WEIGht %.3f",varlist$NOE_DT[[stage]][i,weight])
      NOE_string[idx+1] = sprintf("ASSIgn (resid %i and name %s)(resid %i and name %s) %.1f %.1f %.1f",
                                  varlist$NOE_DT[[stage]][i,Pos1],varlist$NOE_DT[[stage]][i,atom1],
                                  varlist$NOE_DT[[stage]][i,Pos2],varlist$NOE_DT[[stage]][i,atom2],
                                  varlist$NOE_DT[[stage]][i,dist],varlist$NOE_DT[[stage]][i,lower_dist],varlist$NOE_DT[[stage]][i,upper_dist])
      if (varlist$NOE_DT[[stage]][i,!is.na(Pos1_opt2)]) {
        NOE_string[idx+2] = sprintf("OR (resid %i and name %s)(resid %i and name %s)",
                                    varlist$NOE_DT[[stage]][i,Pos1_opt2],varlist$NOE_DT[[stage]][i,atom1_opt2],
                                    varlist$NOE_DT[[stage]][i,Pos2_opt2],varlist$NOE_DT[[stage]][i,atom2_opt2])
        idx=idx+3
      } else {
        idx=idx+2
      }
      varlist$NOE_DT[[stage]][i,restraint_NR := restraint_idx]
      restraint_idx = restraint_idx+1
    }
  }
  
  NOE_string[idx] = "set echo=true end"
  NOE_string[idx+1] = "set wrnlev=5 end"
  
  write(x = NOE_string,file = paste0(varlist$folder,varlist$contacts_noe_file,"_stage",stage,".tbl"))
  return(varlist)
}

write_SS_dihe_restraints = function(varlist,stage) {
  
  #this file writes the dihedral angle restraints derived from secondary structure elements
  
  dihe_string = "set echo=false end"
  dihe_string[2] = "set wrnlev=0 end"
  idx=3
  
  ### set phi
  if (varlist$DIHE_DT[[stage]][type=="phi" & weight > 0.33,.N > 0]) {
    dihe_string[idx] = "!phi"
    idx=idx+1
    for (i in varlist$DIHE_DT[[stage]][,which(type=="phi" & weight > 0.1)]) {
      # dihe_string[idx] = sprintf("WEIGht %.3f",varlist$DIHE_DT[[stage]][i,weight])
      dihe_string[idx] = sprintf("ASSIgn (resid %i and name c) (resid %i and name n) (resid %i and name ca) (resid %i and name c) %.1f %.1f %.1f %i",
                                 varlist$DIHE_DT[[stage]][i,position]-1,varlist$DIHE_DT[[stage]][i,position],
                                 varlist$DIHE_DT[[stage]][i,position],varlist$DIHE_DT[[stage]][i,position],
                                 1,varlist$DIHE_DT[[stage]][i,angle],varlist$DIHE_DT[[stage]][i,delta_angle],2)
      idx = idx+1
    }
  }
  
  ### set psi
  if (varlist$DIHE_DT[[stage]][type=="psi" & weight > 0.33,.N > 0]) {
    dihe_string[idx] = "!psi"
    idx=idx+1
    for (i in varlist$DIHE_DT[[stage]][,which(type=="psi" & weight > 0.1)]) {
      # dihe_string[idx] = sprintf("WEIGht %.3f",varlist$DIHE_DT[[stage]][i,weight])
      dihe_string[idx] = sprintf("ASSIgn (resid %i and name n) (resid %i and name ca) (resid %i and name c) (resid %i and name n) %.1f %.1f %.1f %i",
                                 varlist$DIHE_DT[[stage]][i,position],varlist$DIHE_DT[[stage]][i,position],
                                 varlist$DIHE_DT[[stage]][i,position],varlist$DIHE_DT[[stage]][i,position]+1,
                                 1,varlist$DIHE_DT[[stage]][i,angle],varlist$DIHE_DT[[stage]][i,delta_angle],2)
      idx = idx+1
    }
  }
  
  dihe_string[idx] = "set echo=true end"
  dihe_string[idx+1] = "set wrnlev=5 end"
  
  write(x = dihe_string,file = paste0(varlist$folder,varlist$ss_dihe_file,"_stage",stage,".tbl"))
}

modify_XPLOR_phython = function(varlist,stage) {
  #modify the template phython scripts for XPLOR
  if (stage < 3) {
    #read script
    script_py = scan(file=paste0(varlist$cluster_dir,varlist$protein,"/anneal_template.py"),what="character",sep="\n")
    
    script_py[29] = sprintf("seqToPSF('%s')",varlist$protein_seq)              
    script_py[67] = sprintf("for (name,scale,file) in [('noe_contacts',1,\"%s\"),",paste0(varlist$folder,varlist$contacts_noe_file,"_stage",stage,".tbl"))
    if (varlist$NOE_pot_soft[stage] == T) {
      script_py[72] = "    pot.setPotType(\"soft\") # if you think there may be bad NOEs "    
    }
    script_py[73] = "" #get rid of AveType shortest, default is average
    script_py[85] = sprintf("dihedralRestraintFilename=\"%s\"",paste0(varlist$folder,varlist$ss_dihe_file,"_stage",stage,".tbl"))
    script_py[232] = sprintf("              averageTopFraction=%f,",varlist$top_avg_fraction[stage])
    script_py[233] = "              "
    script_py[236] = sprintf("              averageFilename=\"%sXPLOR_output/SCRIPT_ave.pdb\",",varlist$folder)
    
  } else if (stage == 3) {
    script_py = scan(file=paste0(varlist$cluster_dir,varlist$protein,"/refine_template.py"),what="character",sep="\n")
    script_py[28] = sprintf("protocol.loadPDB(\"%sXPLOR_output/%s\",deleteUnknownAtoms=True)",varlist$folder,paste0(varlist$filename[2],"_ave.pdb"))
    script_py[53] = sprintf("for (name,scale,file) in [('noe_contacts',1,\"%s\"),",paste0(varlist$folder,varlist$contacts_noe_file,"_stage",stage,".tbl"))
    if (varlist$NOE_pot_soft[stage] == T) {
      script_py[57] = "    pot.setPotType(\"soft\") # if you think there may be bad NOEs "    
    }
    script_py[58] = "" #get ride of AveType shortest, default is average
    script_py[69] = sprintf("protocol.initDihedrals(\"%s\"",paste0(varlist$folder,varlist$ss_dihe_file,"_stage",stage,".tbl"))
    
    script_py[82] = sprintf("                    \"resid 1:%d\") # selection should exclude disordered tails",varlist$protein_length)
    
    script_py[228] = sprintf("              averageTopFraction=%f,",varlist$top_avg_fraction[stage]) 
    script_py[231] = sprintf("              averageFilename=\"%sXPLOR_output/SCRIPT_ave.pdb\",",varlist$folder)
  }
  
  script_py[19] = sprintf("outFilename = \"%sXPLOR_output/SCRIPT_STRUCTURE.pdb\"",varlist$folder)
  script_py[20] = sprintf("numberOfStructures=%i",ceiling(varlist$numberOfStructures[stage]))
  
  write(script_py,file=paste0(varlist$folder,varlist$filename[stage],".py"))
}

read_PDBfiles_violations= function(varlist,stage) {
  #read all pdb filenames from current stage
  anneal_pdb = dir(path = paste0(varlist$folder,"XPLOR_output"),pattern = paste0("^",varlist$filename[stage],".*[0-9]*\\.pdb$"))
  for (i in seq_along(anneal_pdb)) {
    pdb_text = scan(paste0(varlist$folder,"XPLOR_output/",anneal_pdb[i],collapse=""), what="character", sep="\n",quiet=T)
    #find lines that give energy term summaries
    pdb_text = pdb_text[grep(pdb_text,pattern="summary")]
    #remove whitspace
    pdb_text = gsub(pattern = "\\s+",replacement = " ",x = pdb_text)
    helper = fread(input = paste0(pdb_text,collapse=" \n "),header=T,fill=T)
    eval(parse(text=paste0("helper_DT = data.table(",paste0(helper[1:.N-1,Name],"=",helper[1:.N-1,Energy],collapse=", "),")")))
    helper_DT[,refRMSD := helper[.N,RMS]]
    helper_DT[,pdb_id := anneal_pdb[i]]
    
    #read TM score and TM RMSD calc
    TM_file=system(command = paste0("TMscore ",varlist$pdb_file," ",varlist$folder,"XPLOR_output/",anneal_pdb[i]),inter=T)
    helper_DT[,TM_score := as.numeric(unlist(strsplit(TM_file[grep(TM_file,pattern="^TM-score")],split = "\\s+"))[3])]
    helper_DT[,TM_RMSD := as.numeric(unlist(strsplit(TM_file[grep(TM_file,pattern="^RMSD")],split = "\\s+"))[6])]
    
    if (i == 1) {
      varlist$energy_XPLOR[[stage]] = helper_DT
    } else {
      varlist$energy_XPLOR[[stage]] = rbind(varlist$energy_XPLOR[[stage]],helper_DT)
    }
  }
  
  
  #read violation files
  anneal_viols = dir(path = paste0(varlist$folder,"XPLOR_output/"),pattern = paste0("^",varlist$filename[stage],".*[0-9]*\\.pdb.viols$"))
  viols_matrix = c()
  varlist$violations[[stage]] = data.table(struct_NR = as.numeric(), NOE_DIHE=as.character(), restraint_NR = as.numeric(), 
                                           Pos1 = as.numeric(), Pos2 = as.numeric(), dev = as.numeric(), type = as.character(), pdb_id = as.character())
  for (i in seq_along(anneal_viols)) {
    viols_text = scan(paste0(varlist$folder,"XPLOR_output/",anneal_viols[i],collapse=""), what="character", sep="\n",quiet=T)
    
    ### extract NOE violations
    #find line-range that gives noe distance violations
    idx_start = grep("Violated NOE restraints in potential term: noe_contacts",viols_text) + 3
    idx_end = grep("number of restraints:",viols_text[idx_start:length(viols_text)]) + idx_start - 2
    #extract lines that give restraint numbers and violation values
    viols = grep("[0-9]+\\s\\(",viols_text[idx_start:idx_end[1]])
    for (j in seq_along(viols)) {
      helper = trimws(viols_text[idx_start+viols[j]-1])
      helper = unlist(strsplit(gsub(pattern = "[()]",replacement = "",x = helper),"\\s+"))
      
      viols_matrix = rbind(viols_matrix,c(i,as.numeric(helper[1]),as.numeric(helper[2]),as.numeric(helper[5]),as.numeric(helper[10])))
      varlist$violations[[stage]] = rbind(varlist$violations[[stage]],data.table(struct_NR = i,NOE_DIHE="NOE",restraint_NR = as.numeric(helper[1])+1,
                                                                                 Pos1 = as.numeric(helper[2]), Pos2 = as.numeric(helper[5]),
                                                                                 dev = as.numeric(helper[10]),type = "",pdb_id = ""))
    }
    
    ### extract NOE violations
    idx_start = grep("   Violations in XPLOR term     CDIH",viols_text)
    idx_ends = grep(" ========================================",viols_text)
    idx_ends = idx_ends[idx_ends > idx_start]
    if (length(idx_ends)>0) {
      for(j in seq_along(idx_ends)) {
        helper = trimws(viols_text[(idx_ends[j]+1):(idx_ends[j]+5)])
        helper2 = c()
        for (k in 1:2) {
          helper2[k] = as.integer(paste0(strsplit(helper[k],split="")[[1]][1:3],collapse=""))
        }
        if (helper2[1] < helper2[2]) {angle_type = "phi"} else {angle_type = "psi"}
        helper3=strsplit(helper[5],split="")[[1]]
        dev = as.numeric(paste0(helper3[(gregexpr(helper[5],pattern = "Delta")[[1]][1]+6) : length(helper3)],collapse=""))
        if (length(varlist$DIHE_DT[[stage]][,which(position == helper2[2] & type == angle_type)])>0) { #sometimes XPLOR seems to generate weird restraints by itself?!
          varlist$violations[[stage]] = rbind(varlist$violations[[stage]],data.table(struct_NR = i,NOE_DIHE="DIHE",
                                                                                     restraint_NR = varlist$DIHE_DT[[stage]][,which(position == helper2[2] & type == angle_type)],
                                                                                     Pos1 = helper2[2], Pos2 = helper2[2],
                                                                                     dev = dev/5,type = angle_type, pdb_id = ""))
        }
      }
    }
  }  
  varlist$violations[[stage]][,pdb_id := anneal_pdb[struct_NR],struct_NR]
  
  ####################
  # #cluster structures by NOE violations and find cluster with lowest total energy
  if (nrow(varlist$violations[[stage]]) > 0) { #if there are any violations....
    res_matrix = matrix(0,nrow=length(anneal_viols),ncol=max(varlist$violations[[stage]]$restraint_NR))
    res_matrix[cbind(varlist$violations[[stage]]$struct_NR,varlist$violations[[stage]]$restraint_NR)] = varlist$violations[[stage]]$dev
    require(pheatmap)
    pdf(paste0(varlist$folder,"R_output/violation_heatmap_stage",stage,".pdf"))
    dev_heatmap = pheatmap(res_matrix,kmeans_k = 4,cluster_cols=F)
    dev.off()
    
    varlist$energy_XPLOR[[stage]][,cluster:= as.factor(dev_heatmap$kmeans$cluster)]
    varlist$energy_XPLOR[[stage]][.N,cluster:= as.factor("avg")]
    
    #find clusters with minimal total energy (use as many clusters as necessary to get to 10% of structures)
    varlist$energy_XPLOR[[stage]][,':=' (cluster_rank = rank(total),N=.N),cluster]
    min_cluster = varlist$energy_XPLOR[[stage]][cluster_rank <= ceiling(varlist$numberOfStructures[stage] * varlist$top_avg_fraction[stage]) & cluster != "avg",
                                                .(mean_total = mean(total),N=unique(N)),cluster][order(mean_total)][min(which(cumsum(N) >= ceiling(varlist$numberOfStructures[stage] * varlist$top_avg_fraction[stage]))),cluster]
    
    min_cluster_frac = ceiling(varlist$numberOfStructures[stage] * varlist$top_avg_fraction[stage])/varlist$energy_XPLOR[[stage]][cluster %in% min_cluster,.N]
    varlist$energy_XPLOR[[stage]][,min_cluster_top10 := FALSE]
    varlist$energy_XPLOR[[stage]][cluster %in% min_cluster & total <=  quantile(varlist$energy_XPLOR[[stage]][cluster %in% min_cluster,total],min_cluster_frac),min_cluster_top10 := TRUE]
    
    #define frequency of restraint violation in the top structures
    top_cluster_viol = unique(varlist$violations[[stage]][pdb_id %in% varlist$energy_XPLOR[[stage]][min_cluster_top10==T,pdb_id],
                                                          .(frac_viol_cluster=.N/varlist$energy_XPLOR[[stage]][min_cluster_top10==T,.N],
                                                            Pos1,Pos2,type),.(restraint_NR,NOE_DIHE)])
    
    # write these to the NOE and DIHE data.tables for the next stage
    if (top_cluster_viol[NOE_DIHE=="NOE",.N>0]) {
      varlist$NOE_DT[[stage+1]] = merge(varlist$NOE_DT[[stage]][,.(weight,Pos1,Pos2,atom1,atom2,dist,lower_dist,upper_dist,
                                                                   Pos1_opt2,Pos2_opt2,atom1_opt2,atom2_opt2,type,restraint_NR)],
                                        top_cluster_viol[NOE_DIHE=="NOE",.(Pos1,Pos2,frac_viol_cluster,restraint_NR)],by=c("Pos1","Pos2","restraint_NR"),all=T)
      varlist$NOE_DT[[stage+1]][is.na(frac_viol_cluster),frac_viol_cluster:=0]
    } else {
      varlist$NOE_DT[[stage+1]] = varlist$NOE_DT[[stage]]
      varlist$NOE_DT[[stage+1]] = cbind(varlist$NOE_DT[[stage+1]],data.table(frac_viol_cluster = 0))
    }
    if (top_cluster_viol[NOE_DIHE=="DIHE",.N>0]) {
      varlist$DIHE_DT[[stage+1]] = merge(varlist$DIHE_DT[[stage]][,.(weight,position,angle,delta_angle,type,ss)],
                                         top_cluster_viol[NOE_DIHE=="DIHE",.(position=Pos1,type, frac_viol_cluster)],by=c("position","type"),all=T)
      varlist$DIHE_DT[[stage+1]][is.na(frac_viol_cluster),frac_viol_cluster:=0]
    } else {
      varlist$DIHE_DT[[stage+1]] = varlist$DIHE_DT[[stage]]
      varlist$DIHE_DT[[stage+1]] = cbind(varlist$DIHE_DT[[stage+1]],data.table(frac_viol_cluster = 0))
    }
  } else {
    varlist$NOE_DT[[stage+1]] = varlist$NOE_DT[[stage]]
    varlist$NOE_DT[[stage+1]] = cbind(varlist$NOE_DT[[stage+1]],data.table(frac_viol_cluster = 0))
    varlist$DIHE_DT[[stage+1]] = varlist$DIHE_DT[[stage]]
    varlist$DIHE_DT[[stage+1]] = cbind(varlist$DIHE_DT[[stage+1]],data.table(frac_viol_cluster = 0))
  }
  ####################
  
  ############ reweight distance constraints for next stage ############ 
  if (stage < 3) {
    varlist$NOE_DT[[stage+1]] = cbind(varlist$NOE_DT[[stage+1]],data.table(weight0 = 0))
    varlist$NOE_DT[[stage+1]][,weight0 := weight]
    varlist$NOE_DT[[stage+1]][,weight := weight0 * (1-frac_viol_cluster)^2]
    
    varlist$DIHE_DT[[stage+1]] = cbind(varlist$DIHE_DT[[stage+1]],data.table(weight0 = 0))
    varlist$DIHE_DT[[stage+1]][,weight0 := weight]
    varlist$DIHE_DT[[stage+1]][,weight := weight0 * (1-frac_viol_cluster)^2]
  }
  
  return(varlist)
} 

