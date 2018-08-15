
######################################################################################
### extract position-pair distances and secondary structure elements from PDB file ###
######################################################################################

pairdistances_from_PDB = function(input_file,
                                  dataset_dir,
                                  given_chainids = "A",
                                  aa_seq,
                                  idx_pdb_start = 1,
                                  idx_DMS_start = 1,
                                  idx_DMS_end = NA,
                                  dist_cutoff = 8,
                                  debug_this=F,
                                  suffix = "") {
  
  ### variables 
  # input_file: PDB file
  # dataset_dir: dataset directory, like "GB1/", 
  #### > will deposit two .txt file (on for pair distances and one for secondary structure) to dataset_dir/processed_data/ and plot contactmaps to /dataset_dir/results/preprocessing/
  # given_chainids: which chain in the PDB file to extract distances from; in case of protein-protein interactions needs to be a vector with both chains that trans-distances should be calculated over, e.g. c("A","C")
  # aa_seq: amino acid sequence of reference structure (in DMS data)
  # idx_pdb_start: first position in PDB to consider
  # idx_DMS_start: first position in reference sequence to consider
  # idx_DMS_end: first position in reference sequence to consider, if NA (default), it will compare the full reference sequence to the PDB file sequence
  # dist_cutoff: in Angstrom, used for plotting the contact map
  # debug_this: if TRUE, the function will stop after printing comparision between PDB seq and DMS seq to adjust position indicies if necessary
  # suffix: to be added to 
  
  
  require(data.table)
  require(Rpdb) 
  require(pdist)
  require(ggplot2)
  require(cowplot)
  
  #if idx_DMS_end is not given, compare across full sequence length
  if (is.na(idx_DMS_end)) {
    idx_DMS_end = c(0)
    for (c in seq_along(given_chainids)) {
      idx_DMS_end[c] = length(strsplit(aa_seq[[c]],"")[[1]])
    }
  }
  
  #load PDB structure
  PDB_structure = read.pdb(input_file,MODEL=NULL)
  
  #for PDB files with NMR ensembles evaluate each model
  M = length(grep("MODEL",names(PDB_structure)))
  if (M==0) {M=1}
  for (m in 1:M) {
    
    #load model to structure data.table
    if (length(grep("MODEL",names(PDB_structure)))==0) {
      structure = data.table(eval(parse(text=paste0("PDB_structure$atoms"))))
    } else {
      structure = data.table(eval(parse(text=paste0("PDB_structure$MODEL.",m,"$atoms"))))
    }
    
    #restrict to ATOM entries
    structure = structure[recname=="ATOM" & chainid %in% given_chainids]
    
    #extract amino acid sequence from PDB file
    aaseq_PDB = unique(structure[,.(AA = convert_AAabr_one_three(as.character(unique(resname))),chainid),by=resid])
    setkey(aaseq_PDB,resid)
    
    ## compare given DMS aaseq and aaseq from PDB file
    DMS_aa_seq = list()
    PDB_aa_seq = list()
    for (c in seq_along(given_chainids)) {
      DMS_aa_seq[[c]] = strsplit(aa_seq[[c]],"")[[1]][idx_DMS_start[c]:idx_DMS_end[c]]
      PDB_aa_seq[[c]] = aaseq_PDB[chainid == given_chainids[c]][.(idx_pdb_start[c]:(idx_pdb_start[c]+length(DMS_aa_seq[[c]])-1)),AA]
      PDB_aa_seq[[c]][is.na(PDB_aa_seq)] = "X"
      
      if (m==1) {
        print(paste0(input_file," chain ",given_chainids[c]))
        print(paste0('DMS seq [',idx_DMS_start[c],':',idx_DMS_end[c],'] ',paste0(DMS_aa_seq[[c]],collapse="")))
        print(paste0(sum(DMS_aa_seq[[c]]==PDB_aa_seq[[c]]),'/',length(DMS_aa_seq[[c]]),'          ',paste0(as.numeric(DMS_aa_seq[[c]] == PDB_aa_seq[[c]]),collapse = "")))
        print(paste0('PDB seq [',idx_pdb_start[c],':',(idx_pdb_start[c]+length(DMS_aa_seq[[c]])-1),'] ',paste0(PDB_aa_seq[[c]],collapse="")))
      }
    }
    
    #if debug_this ==T function will stop here to make adjustments to position indicies
    if (debug_this) {
      browser()
    }
    
    #initialize indicies for distance calcualtions
    if (length(given_chainids) == 1) {
      two_chainids = rep(given_chainids,2)
      two_starts_DMS = rep(idx_DMS_start,2)
      two_starts = rep(idx_pdb_start,2)
      two_ends = rep(idx_pdb_start + (idx_DMS_end-idx_DMS_start),2)
    } else {
      two_chainids = given_chainids
      two_starts_DMS = idx_DMS_start
      two_starts = idx_pdb_start
      two_ends = idx_pdb_start + (idx_DMS_end-idx_DMS_start)
    }
    
    #initialize distance table
    if (m==1) {
      distance = data.table(Pos1=rep(two_starts[1]:two_ends[1],two_ends[2]-two_starts[2]+1),
                            Pos2=rep(two_starts[2]:two_ends[2],each=two_ends[1]-two_starts[1]+1))
      setkey(distance,Pos1,Pos2)
      distance[,WT_AA1 := convert_AAabr_one_three(as.character(unique(structure[chainid == two_chainids[1] & resid == Pos1,resname]))),Pos1]
      distance[,WT_AA2 := convert_AAabr_one_three(as.character(unique(structure[chainid == two_chainids[2] & resid == Pos2,resname]))),Pos2]
      distance[,chainids := paste0(given_chainids)]
    }
    
    #calculate minimal side-chain heavy atom distance
    structure_HA = structure[union(intersect(grep(pattern="^[COSN][B-Z]$",elename),
                                             which(resname != "GLY")),intersect(grep(pattern="^CA$",elename),which(resname == "GLY")))]
    distance[,paste0("scHAmin",m):=min(as.matrix(pdist(as.matrix(structure_HA[chainid == two_chainids[1] & resid == Pos1,.(x1,x2,x3)]),
                                                       as.matrix(structure_HA[chainid == two_chainids[2] & resid == Pos2,.(x1,x2,x3)])))),
             by=.(Pos1,Pos2)]
    
    #calculate minimal all heavy atom distance
    structure_HA = structure[!grepl(pattern="H",elename)]
    distance[,paste0("HAmin",m):=min(as.matrix(pdist(as.matrix(structure_HA[chainid == two_chainids[1] & resid == Pos1,.(x1,x2,x3)]),
                                                     as.matrix(structure_HA[chainid == two_chainids[2] & resid == Pos2,.(x1,x2,x3)])))),
             by=.(Pos1,Pos2)]
    
    #calculate CB distances (use CA in case of Glycine)
    structure_CB = structure[elename == "CB" | (elename == "CA" & resname == "GLY"),.(chainid,resid,x1,x2,x3)]
    distance[,paste0("CB",m):=min(as.matrix(pdist(as.matrix(structure_CB[chainid == two_chainids[1] & resid == Pos1,.(x1,x2,x3)]),
                                                  as.matrix(structure_CB[chainid == two_chainids[2] & resid == Pos2,.(x1,x2,x3)])))),
             by=.(Pos1,Pos2)]
  }
  
  ## average over distances
  distance[,scHAmin := rowMeans(.SD),by=.(Pos1,Pos2),.SDcols = grep("scHAmin[0-9]",names(distance))]
  distance[,HAmin := rowMeans(.SD),by=.(Pos1,Pos2),.SDcols = grep("HAmin[0-9]",names(distance))]
  distance[,CB := rowMeans(.SD),by=.(Pos1,Pos2),.SDcols = grep("CB[0-9]",names(distance))]
  
  
  #if there's multiple structural models, average over all and also calculate uncertainity
  if (length(grep("MODEL",names(PDB_structure)))==0) {
    distance[,scHAmin_sd := 0]
    distance[,HAmin_sd := 0]
    distance[,CB_sd := 0]
  } else {
    distance[,scHAmin_sd := stats::sd(.SD),by=.(Pos1,Pos2),.SDcols = grep("scHAmin[0-9]",names(distance))]
    distance[,HAmin_sd := stats::sd(.SD),by=.(Pos1,Pos2),.SDcols = grep("HAmin[0-9]",names(distance))]
    distance[,CB_sd := stats::sd(.SD),by=.(Pos1,Pos2),.SDcols = grep("CB[0-9]",names(distance))]
  } 
  
  ##adjust positions to positions in alignment DMS/PDB sites
  contactmap = distance[between(Pos1,two_starts[1],two_ends[1]) & 
                          between(Pos2,two_starts[1],two_ends[2]),
                        .(Pos1 = Pos1 - (two_starts[1] - two_starts_DMS[1]),
                          Pos2 = Pos2 - (two_starts[2] - two_starts_DMS[2]),
                          WT_AA1,WT_AA2,chainids,
                          HAmin,scHAmin,CB,HAmin_sd,scHAmin_sd,CB_sd)]
  
  #save pairwise distance table
  pdb_filename = strsplit(strsplit(input_file,"/")[[1]][length(strsplit(input_file,"/")[[1]])],"\\.")[[1]][1]
  write.table(file = paste0(dataset_dir,'processed_data/PDB_contactmap_',pdb_filename,'_',paste0(given_chainids,collapse=""),suffix,".txt",collapse = ""),
              x = contactmap,quote = F,row.names = F,col.names = T)
  
  
  ################################################# 
  ### extract secondary structure from PDB file ###
  #################################################
  if (length(given_chainids) == 1) { #if looking at a single chain
    secondary_structure = data.table(Pos = idx_DMS_start:idx_DMS_end,ss = "C")
    
    output = scan(file=input_file,what="character",sep="\n")
    
    helix = output[grep(output,pattern="^HELIX")]
    if (length(helix)>0) {
      helix1 = sapply(X=1:length(helix),FUN = function(X){strsplit(helix[X],split="\\s+")[[1]]})
      helix2 = data.table(t(helix1[4:9,]))
      names(helix2) = c("aa1","chainid1","pos1","aa2","chainid2","pos2")
      helix3 = helix2[chainid1 == given_chainids]
      if (nrow(helix3)>0) {
        for (i in 1:nrow(helix3)) {
          secondary_structure[between(Pos,
                                      helix3[i,as.numeric(pos1) - (two_starts[1] - two_starts_DMS[1])],
                                      helix3[i,as.numeric(pos2) - (two_starts[1] - two_starts_DMS[1])]),
                              ss := "H"]
        }
      }
    }
    
    strand = output[grep(output,pattern="^SHEET")]
    if (length(strand)>0) {
      strand1 = sapply(X=1:length(strand),FUN = function(X){strsplit(strand[X],split="\\s+")[[1]]})  
      if (is.list(strand1)) {
        strand2 = data.table(t(strand1[[1]][5:10]))
        for (l in 2:length(strand1)) {
          strand2 = rbind(strand2,data.table(t(strand1[[l]][5:10])))
        }
      } else {
        strand2 = data.table(t(strand1[5:10,]))  
      }
      
      names(strand2) = c("aa1","chainid1","pos1","aa2","chainid2","pos2")
      strand3 = strand2[chainid1 == given_chainids]
      if (nrow(strand3) > 0) {
        for (i in 1:nrow(strand3)) {
          secondary_structure[between(Pos,
                                      strand3[i,as.numeric(pos1) - (two_starts[1] - two_starts_DMS[1])],
                                      strand3[i,as.numeric(pos2) - (two_starts[1] - two_starts_DMS[1])]),
                              ss := "E"]
        }
      }
    }
    names(secondary_structure)[2] = "PDB"
    
    write.table(paste0(dataset_dir,"processed_data/PDB_secondary_structure_",pdb_filename,"_",given_chainids,suffix,".txt"),
                x = secondary_structure,quote = F,row.names = F,col.names = T)
    
    #for plotting
    secondary_structure[,rleidx := rleid(PDB)]
  }
  
  
  
  theme_set(theme_classic(base_size=9))
  #plot contact map
  P1=ggplot() +   
    geom_raster(data=contactmap,aes(x=Pos1,y=Pos2,fill=HAmin<dist_cutoff),show.legend = F) +
    scale_fill_manual(values = c("white","grey")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_reverse(expand = c(0,0)) +
    labs(fill = "<HAmin>",title=paste("all heavy atom distance < ",dist_cutoff,"A"),x="Pos1",y="Pos2")
  if (length(given_chainids) == 1) { #add secondary structure
    P1 = P1 +  
      geom_segment(data = secondary_structure[,.(start = min(Pos)-0.5,end = max(Pos)+0.5,ss=unique(PDB)),rleidx],
                   aes(x=start,y=start,xend=end,yend=end,color=ss,size=ss),show.legend = F) +
      scale_size_manual(breaks = c("C","H","E"),values = c(0.5,1.5,1.5)) +
      scale_color_manual(breaks = c("C","H","E"),values = c("black","orange","darkgreen"))
  }
  
  P2=ggplot() +   
    geom_raster(data=contactmap,aes(x=Pos1,y=Pos2,fill=scHAmin<dist_cutoff),show.legend = F) +
    scale_fill_manual(values = c("white","grey")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_reverse(expand = c(0,0)) +
    labs(fill = "<scHAmin>",title=paste("side-chain heavy atom distance < ",dist_cutoff,"A"),x="Pos1",y="Pos2")
  if (length(given_chainids) == 1) { #add secondary structure
    P2 = P2 +  
      geom_segment(data = secondary_structure[,.(start = min(Pos)-0.5,end = max(Pos)+0.5,ss=unique(PDB)),rleidx],
                   aes(x=start,y=start,xend=end,yend=end,color=ss,size=ss),show.legend = F) +
      scale_size_manual(breaks = c("C","H","E"),values = c(0.5,1.5,1.5)) +
      scale_color_manual(breaks = c("C","H","E"),values = c("black","orange","darkgreen"))
  }
  
  P3=ggplot() +   
    geom_raster(data=contactmap,aes(x=Pos1,y=Pos2,fill=CB<dist_cutoff),show.legend = F) +
    scale_fill_manual(values = c("white","grey")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_reverse(expand = c(0,0)) +
    labs(fill = "<CB>",title=paste("Cbeta distances < ",dist_cutoff,"A"),x="Pos1",y="Pos2")
  if (length(given_chainids) == 1) { #add secondary structure
    P3 = P3 +  
      geom_segment(data = secondary_structure[,.(start = min(Pos)-0.5,end = max(Pos)+0.5,ss=unique(PDB)),rleidx],
                   aes(x=start,y=start,xend=end,yend=end,color=ss,size=ss),show.legend = F) +
      scale_size_manual(breaks = c("C","H","E"),values = c(0.5,1.5,1.5)) +
      scale_color_manual(breaks = c("C","H","E"),values = c("black","orange","darkgreen"))
  }
  
  
  #plot this fourth one just for the secondary structure legend
  P4=ggplot() +   
    geom_raster(data=contactmap,aes(x=Pos1,y=Pos2,fill=scHAmin),show.legend = F) +
    scale_fill_gradient(low="white",high="grey") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_reverse(expand = c(0,0)) +
    labs(fill = "",title="scHAmin absolute distances",x="Pos1",y="Pos2")
  if (length(given_chainids) == 1) { #add secondary structure
    P4 = P4 +  
      geom_segment(data = secondary_structure[,.(start = min(Pos)-0.5,end = max(Pos)+0.5,ss=unique(PDB)),rleidx],
                   aes(x=start,y=start,xend=end,yend=end,color=ss,size=ss),show.legend = T) +
      scale_size_manual(breaks = c("C","H","E"),values = c(0.5,1.5,1.5)) +
      scale_color_manual(breaks = c("C","H","E"),values = c("black","orange","darkgreen"))
  }
  
  P = plot_grid(plotlist = list(P1,P2,P3,P4),nrow=2)
  P
  ggsave(plot = P,filename = paste0(dataset_dir,'results/preprocessing/contactmap_',pdb_filename,'_',paste0(given_chainids,collapse=""),suffix,'.pdf',collapse = ""),width = 8,height=7)
}

