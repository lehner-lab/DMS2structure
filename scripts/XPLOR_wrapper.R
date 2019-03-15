
XPLOR_wrapper = function(input_PWI,
                         SS_mode = "SSonly",
                         input_SS_file,
                         L = c(0.5,1,1.5,2),
                         predictorXlength = c(),
                         prefix,
                         dataset_dir,
                         protein_sequence,
                         pdb_file,
                         cores = 12,
                         queue = "long-sl7,short-sl7",
                         numberOfStructures = c(500,500,500),
                         top_avg_fraction = c(0.1,0.1,0.1),
                         NOE_pot_soft = c(TRUE,TRUE,FALSE),
                         home_dir,
                         cluster_dir,
                         reporting_email = NA,
                         login_serveraddress,
                         debug_this = F,
                         linear_dist = 5,
                         dist_restraint = 8) {
  
  #these scripts are multi-level wrapper functions to run XPLOR-NIH on a compute cluster
  #they locally create/modify the necessary scripts and folder structure, and scp this to the cluster, then remotely execute parallel jobs
  #this might take quite some work to adapt to other working environments
  
  ### variables
  # input_PWI: pairwise interaction score data.table; this should have Pos1, Pos2, WT_AA1 and WT_AA2 columns + all interscores that restraints should be derived from
  # SS_mode: - either "SSonly", in which case it only uses predicted secondary structure elements (from DMS data or e.g. PSIPRED) for restraints
  #          - or "SSsheets", in which it also derives restraints for beta sheet pairing hbonding
  # input_SS_file: secondary structure element input file, 
  #                   - if SS_mode == "SSonly" this is a table with position columns and either ONE ss index column to be used for all given interaction scores, or a SS index column per interaction score (with the interaction score as column name)
  #                   - if SS_mode == "SSsheets" this needs to be a secondary structure element RData file that also contains the beta sheet pairing table
  # L: number of top contacts * protein length used for tertiary contact restraints, can be a vector, a simulation per score x L combination will be started [additionally, if there is a "control" score, the script will create a negative control, i.e. L = 0]
  # predictorXlength: this can be a data.table with first column indicating interaction score, second column indicate L; this will rerun only the specific score:L combinations indicated instead of all scores versus all L + negative control
  # dataset_dir: dataset directory, like "GB1/", it will put results dataset_dir/results/PWI/
  # prefix: to be added to results files (in case of running diff. versions of data from same dataset etc)
  # protein_sequence: amino acid sequence of the protein to be modeled, used by XPLOR
  # pdb_file (optional): reference structure pdb file to compare structural models to in terms of RMSD and template modeling score
  # cores: number of cores to request on the computing cluster, max is 16
  # queue: which queues to submit jobs to
  # numberOfStructures: how many structures to create in each of the three modeling stages
  # top_avg_fraction: fraction of top evaluated models (by XPLOR total energy) to use to decide restraint violations (stages 1+2) and calculate an average structural model (stage 2 as starting point for stage 3, and stage 3 as a final output)
  # NOE_pot_soft: for the three different stages; TRUE: use soft well potentials for distance restraints (if potentially many restraitns from false positive contacts)
  # home_dir: a directory on local machine in which a "tmp" folder is create to build the script/folder structure locally then copy to cluster; this is removed at the end of the script; could be "/Users/me/"
  # cluster_dir: (main) directory on the compute cluster into which the folder structure should be copied and where the simulations are executed
  # reporting_email: if not NA: email to report to about job status for qsub system
  # login_serveraddress: mylogin@serveraddress.com
  # debug_this: if TRUE, script will stop at certain points for debugging
  # linear_dist: minimum linear (residue) distance for defining top predicted contacts
  # dist_restraint: distance restraint (in Angstrom) for top predicted contacts
  
  require(data.table)
  require(ssh.utils)
  #load utility scripts
  filelist = list.files('scripts/')
  sapply(paste0('scripts/',filelist),source,.GlobalEnv)
  
  #create tmp directory for data structure
  system(command = paste0("mkdir ",home_dir,"tmp"),wait=T)
  system(command = paste0("mkdir ",home_dir,"tmp/",prefix),wait=T)
  
  #which interaction score to use as predictors?
  predictor = setdiff(names(input_PWI),c("Pos1","Pos2","WT_AA1","WT_AA2"))
  #create data.table with predictor x L combinations
  if (length(predictorXlength) == 0) {
    predictorXlength = data.table(expand.grid(predictor,L))
    names(predictorXlength) = c("predictor","L")
    if (length(grep("control",predictor)) > 0) { #if preditors include controls, also create a negative control without tertiary distance restraints
      predictorXlength = rbind(predictorXlength,data.table(predictor = "control",L = 0))
    }
  } else { #the table was supplied with specific combinations
    names(predictorXlength) = c("predictor","L")
  }
  predictorXlength[,predictor := as.character(predictor)]
  
  if (debug_this) {browser()}
  
  #set up list structure for variables and results
  varlist = list()
  varlist$cores = cores
  varlist$protein = prefix
  varlist$cluster_dir = cluster_dir
  varlist$numberOfStructures = numberOfStructures
  varlist$top_avg_fraction = top_avg_fraction
  varlist$NOE_pot_soft = NOE_pot_soft
  varlist$filename = c("anneal_stage1","anneal_stage2","refine")
  varlist$contacts_noe_file = "NOE_restraints"
  varlist$ss_dihe_file = "DIHE_restraints"
  
  
  varlist$protein_length = nchar(protein_sequence)
  ### define .seq file for XPLOR, 3letter AA names 
  write(x = paste(toupper(sapply(strsplit(protein_sequence,"")[[1]],convert_AAabr_one_three)),collapse=" "), 
        file=paste0(home_dir,"tmp/",prefix,"/protein.seq"))
  varlist$protein_seq = paste0(cluster_dir,prefix,"/protein.seq") #replace by fasta seq, then create .seq file and copy
  
  #copy pdb template file
  system(command = paste0("cp ",pdb_file," ",home_dir,"tmp/",prefix,"/template.pdb"),wait=T)
  varlist$pdb_file = paste0(cluster_dir,prefix,"/template.pdb")
  
  
  ### per combination of predictor and L, create restraint files etc
  for (idx in 1:nrow(predictorXlength)) {
    varlist$predictor = predictorXlength$predictor[idx]
    varlist$L = predictorXlength$L[idx]
    varlist$folder = paste0(cluster_dir,prefix,"/",varlist$predictor,"/L",varlist$L,"/")
    
    
    ### initiate lists
    #for restraints
    varlist$NOE_DT = list()
    varlist$NOE_DT[[1]] = data.table(weight = as.numeric(),Pos1 = as.numeric(), Pos2 = as.numeric(), 
                                     atom1 = as.character(), atom2 = as.character(), 
                                     dist = as.numeric(), lower_dist = as.numeric(), upper_dist = as.numeric(),
                                     Pos1_opt2 = as.numeric(), Pos2_opt2 = as.numeric(),
                                     atom1_opt2 = as.character(),atom2_opt2 = as.character(),type = as.character())
    varlist$DIHE_DT = list()
    varlist$DIHE_DT[[1]] = data.table(weight = as.numeric(),position = as.numeric(),angle=as.numeric(), 
                                      delta_angle=as.numeric(),type = as.character(),ss = as.character())
    #for evaluation
    varlist$energy_XPLOR = list()
    varlist$violations = list()
    
    ##### secondary structure elements
    if (SS_mode == "SSonly") {
      ### only secondary structure
      if (varlist$predictor == "control") {
        files = list.files(path = paste0(dataset_dir,"processed_data/"))
        load(paste0(dataset_dir,"processed_data/",files[grep("secondary_structure_elements_control.*RData",files)]))
        input_SS = ss_elements
      } else {
        input_SS = fread(paste0(dataset_dir,"processed_data/",input_SS_file))
      }
      ### define secondary structure restraints
      if (ncol(input_SS) > 2){
        varlist$SS_pred = input_SS[,.(position = Pos,ss = .SD),,.SDcols = varlist$predictor]
      } else if (ncol(input_SS) == 2) {
        varlist$SS_pred = input_SS[,.(position = Pos,ss = .SD),,.SDcols = 2]
      } else {
        print("number secondary structure inputs doesn't match the number of input features")
        secondary_structure = error
      }
      varlist$SS_pred[,rleidx := rleid(ss)]
    } else if (SS_mode == "SSsheets") { ## also include beta sheet hbonding
      
      ### secondary structure and beta sheet hbonding
      if (varlist$predictor == "control") {
        files = list.files(path = paste0(dataset_dir,"processed_data/"))
        load(paste0(dataset_dir,"processed_data/",files[grep("secondary_structure_elements_control.*RData",files)]))
        beta_hbonds = hbonds_from_betasheetpairing(beta_sheet_pairing)
      } else {
        load(paste0(dataset_dir,"processed_data/",input_SS_file,"secondary_structure_elements_",varlist$predictor,".RData"))
        beta_hbonds = hbonds_from_betasheetpairing(beta_sheet_pairing,ss_data)
      }
      # secondary structure
      varlist$SS_pred = ss_elements
      names(varlist$SS_pred) = c("position","ss")
      varlist$SS_pred[,rleidx := rleid(ss)]
      # write beta sheet hbonding to restraints list
      if (nrow(beta_hbonds)>0) {
        varlist$NOE_DT[[1]] = rbind(varlist$NOE_DT[[1]], data.table(weight=1, Pos1=beta_hbonds$hn_opt1,Pos2=beta_hbonds$o_opt1, atom1="hn", atom2="o",
                                                                    dist=2,lower_dist=0.2, upper_dist = 0.1,
                                                                    Pos1_opt2=beta_hbonds$hn_opt2,Pos2_opt2=beta_hbonds$o_opt2, 
                                                                    atom1_opt2="hn", atom2_opt2="o",type = "beta sheet hbond"))
        #XPLOR doesn't create an O for last position when using seq2PSF
        varlist$NOE_DT[[1]] = varlist$NOE_DT[[1]][(is.na(Pos2_opt2) & !(Pos2 == varlist$protein_length)) |
                                                    (!is.na(Pos2_opt2) & !(Pos2 == varlist$protein_length | Pos2_opt2 == varlist$protein_length))]
      }
    }
    
    ##### write phi/psi dihedral angle restraints for secondary structure
    # angle values +- delta for PSI and PHI angles in secondary structure elements
    # these values are taken from Table1 of "CONFOLD: Residue-residue contact-guided ab initio protein folding", Adhikari et al. 2015
    #alpha
    phiH=-63.5
    dphiH=4.5
    psiH=-41.5
    dpsiH=5
    #beta
    phiE=-118
    dphiE=10.7
    psiE=134
    dpsiE=8.6
    
    varlist$DIHE_DT[[1]] = rbind(varlist$DIHE_DT[[1]],varlist$SS_pred[ss=="E" & position > 1,
                                                                      .(weight=1,position,angle = phiE,delta_angle = dphiE,type="phi",ss="beta")])
    varlist$DIHE_DT[[1]] = rbind(varlist$DIHE_DT[[1]],varlist$SS_pred[ss=="E" & position < varlist$protein_length,
                                                                      .(weight=1,position,angle = psiE,delta_angle = dpsiE,type="psi",ss="beta")])
    varlist$DIHE_DT[[1]] = rbind(varlist$DIHE_DT[[1]],varlist$SS_pred[ss=="H" & position > 1,
                                                                      .(weight=1,position,angle = phiH,delta_angle = dphiH,type="phi",ss="alpha")])
    varlist$DIHE_DT[[1]] = rbind(varlist$DIHE_DT[[1]],varlist$SS_pred[ss=="H" & position < varlist$protein_length,
                                                                      .(weight=1,position,angle = psiH,delta_angle = dpsiH,type="psi",ss="alpha")])
    
    # if (debug_this) {browser()}
    #write O-O restraints in beta strands
    pos = varlist$SS_pred[,.(ss,next_E = varlist$SS_pred[position == (unlist(.SD)+1),ss=="E"]),
                          position,.SDcols = "position"][ss=="E" & next_E & position < varlist$protein_length-1,position]
    if (length(pos) > 0) {
      varlist$NOE_DT[[1]] = rbind(varlist$NOE_DT[[1]], data.table(weight=1,Pos1=pos,Pos2=pos+1,atom1 = "o", atom2 = "o",
                                                                  dist = 4.5,lower_dist=0.1, upper_dist=0.1,
                                                                  Pos1_opt2 = NA, Pos2_opt2 = NA, atom1_opt2 = NA,atom2_opt2 = NA, type = "beta strand"))
    }
    
    ###### define restraints from top predicted contacts (Cbeta distances) ###### 
    # get scores for predictor (only those > linear_dist apart in linear sequence)
    PWI = input_PWI[Pos1 < Pos2-linear_dist,cbind(Pos1,Pos2,WT_AA1,WT_AA2,.SD),,.SDcols = varlist$predictor]
    names(PWI)[5] = "score"
    
    #exclude contacts within predicted secondary structure
    for (i in varlist$SS_pred[ss %in% c("E","H"),unique(rleidx)]) {
      PWI[Pos1 %in% varlist$SS_pred[rleidx==i,position] & Pos2 %in% varlist$SS_pred[rleidx==i,position],score := NA]
    }
    if (varlist$predictor == "control") {
      set.seed(seed=1603)
      helper = PWI[Pos1<Pos2-linear_dist & score==T][sample(1:.N,((varlist$protein_length-linear_dist)*varlist$L))]
    } else {
      helper = PWI[Pos1<Pos2-linear_dist][order(-score)][1:((varlist$protein_length-linear_dist)*varlist$L)]
    }
    helper[,':=' (atom1 = "cb",atom2 = "cb")]
    helper[WT_AA1 == "G",atom1 := "ca"]
    helper[WT_AA2 == "G",atom2 := "ca"]
    if (nrow(helper) > 0) {
      varlist$NOE_DT[[1]] = rbind(varlist$NOE_DT[[1]], data.table(weight=helper[,score/mean(score)],Pos1 = helper$Pos1,Pos2 = helper$Pos2,
                                                                  atom1 = helper$atom1,atom2 = helper$atom2,
                                                                  dist = dist_restraint, lower_dist=dist_restraint, upper_dist=0,
                                                                  Pos1_opt2 = NA, Pos2_opt2 = NA, atom1_opt2 = NA,atom2_opt2 = NA, type = "contact"))
    }
    
    
    ########## write files
    #create predictor sub-directories
    system(command = paste0("mkdir ",home_dir,"tmp/",prefix,"/",predictorXlength$predictor[idx]),wait=T)
    system(command = paste0("mkdir ",home_dir,"tmp/",prefix,"/",predictorXlength$predictor[idx],"/L",predictorXlength$L[idx]),wait=T)
    #save varlist
    save(file = paste0(home_dir,"tmp/",prefix,"/",varlist$predictor,"/L",varlist$L,"/varlist.RData"),
         list = "varlist")
    
    ##### copy modeling scripts
    system(command = paste0("cp -r scripts/XPLOR/XPLOR_simulations.R ",home_dir,"tmp/",prefix,"/"),wait=T)
    system(command = paste0("cp -r scripts/XPLOR/XPLOR_modeling_functions_v2.R ",home_dir,"tmp/",prefix,"/"),wait=T)
    system(command = paste0("cp -r scripts/XPLOR/anneal_template.py ",home_dir,"tmp/",prefix,"/"),wait=T)
    system(command = paste0("cp -r scripts/XPLOR/refine_template.py ",home_dir,"tmp/",prefix,"/"),wait=T)
    
    
    ######## write bash file to execute qsub job on cluster
    bash_script = "#!/bin/bash"
    bash_script[2] = paste0("#$ -q ",queue)
    bash_script[3] = paste0("#$ -N ",prefix,"_L",varlist$L,"_",varlist$predictor)
    if (!is.na(reporting_email)) {
      bash_script[4] = paste0("#$ -M ",reporting_email)
    } else {bash_script[4] = ""}
    bash_script[5] = "#$ -m ae"
    bash_script[6] = paste0("#$ -l virtual_free=",varlist$cores*3,"G")
    bash_script[7] = paste0("#$ -pe smp ",varlist$cores)
    bash_script[8] = paste0("#$ -o ",varlist$folder)
    bash_script[9] = paste0("#$ -e ",varlist$folder)
    bash_script[10] = paste0("/software/bl/el7.2/R-3.4.0/bin/Rscript --vanilla ",cluster_dir,prefix,"/XPLOR_simulations.R -v ",varlist$folder,"varlist.RData")
    
    write(bash_script,
          file = paste0(home_dir,"tmp/",prefix,"/",varlist$predictor,"/L",varlist$L,"/bash_execute_XPLOR_sim.sh"))
    
  }

###### copy folder structure to server
system(command = paste0("scp -r ",home_dir,"tmp/",prefix," ",login_serveraddress,":",cluster_dir),wait=T)


if (debug_this) {browser()}

###### execute bash_scripts
for (idx in 1:nrow(predictorXlength)) {
  run.remote(paste0("qsub ",cluster_dir,prefix,"/",predictorXlength$predictor[idx],"/L",predictorXlength$L[idx],"/bash_execute_XPLOR_sim.sh"),remote = login_serveraddress)
}


###### remove tmp dir on local machine
system(command = paste0("rm -rf ",home_dir,"tmp/"))
}