
#################################################################
### extract position-pair distances for polymer from PDB file ###
#################################################################

pairdistances_from_PDB_crystal = function(input_file,
  dataset_dir,
  aa_seq,
  idx_pdb_start = 1,
  idx_DMS_start = 1,
  idx_DMS_end = NA,
  dist_cutoff = 8,
  debug_this=F,
  suffix = ""){
  
  ### variables 
  # input_file: PDB file
  # dataset_dir: dataset directory, like "GB1/", 
  #### > will deposit three .txt file (one for pair distances and one for secondary structure of psuedomonomer, one for pair distances of polymer) to dataset_dir/processed_data/ and plot contactmaps to /dataset_dir/results/preprocessing/
  # aa_seq: amino acid sequence of reference (monomer) structure (in DMS data) -- this is compared to inferred monomer sequence from PDB file and produces error if they do not agree
  # idx_pdb_start: first position in PDB to consider (WARNING: argument passed to pairdistances_from_PDB)
  # idx_DMS_start: first position in reference sequence to consider (WARNING: argument passed to pairdistances_from_PDB)
  # idx_DMS_end: first position in reference sequence to consider, if NA (default), it will compare the full reference sequence to the PDB file sequence (WARNING: argument passed to pairdistances_from_PDB)
  # dist_cutoff: in Angstrom, used for plotting the contact map (WARNING: argument passed to pairdistances_from_PDB)
  # debug_this: if TRUE, the function will stop after printing comparision between PDB seq and DMS seq to adjust position indicies if necessary (WARNING: argument passed to pairdistances_from_PDB)
  # suffix: to be added to (WARNING: argument passed to pairdistances_from_PDB)

  require(data.table)

  #Atom lines
  atom_lns <- c("ATOM")
  #Reformat PDB file
  #Read PDB file
  pdb_tab <- read.pdb(input_file)
  #Get all atoms
  temp_atoms <- pdb_tab$atoms[pdb_tab$atoms$recname %in% atom_lns,]
  #Convert structure into "pseudomonomer" (residue ids stricly ascending and same chain)
  resid_rle <- rle(temp_atoms$resid)
  pdb_tab$atoms[pdb_tab$atoms$recname %in% atom_lns,]$resid <- rep(resid_rle$values[1]:(resid_rle$values[1]+length(resid_rle$values)-1), times = resid_rle$lengths)
  pdb_tab$atoms[pdb_tab$atoms$recname %in% atom_lns,]$chainid <- "A"
  #Amino acid sequence of pseudomonomer
  pseudomonomer_length <- length(unique(pdb_tab$atoms[pdb_tab$atoms$recname %in% atom_lns,c("resname", "resid")])$resname)
  aa_seq_pseudomonomer <- paste0(rep(aa_seq, pseudomonomer_length/nchar(aa_seq)), collapse = "")
  #Write to PDB file
  input_file_pseudomonomer <- file.path(dataset_dir, "processed_data", paste0(strsplit(basename(input_file), "\\.")[[1]][1], "_pseudomonomer", suffix, ".pdb"))
  write.pdb(pdb_tab, file = input_file_pseudomonomer)
  #Get pair distances
  pairdistances_from_PDB(input_file_pseudomonomer, dataset_dir = dataset_dir, aa_seq = aa_seq_pseudomonomer, 
    idx_pdb_start = idx_pdb_start, idx_DMS_start = idx_DMS_start, idx_DMS_end = idx_DMS_end, dist_cutoff = dist_cutoff, debug_this = debug_this, suffix = suffix)
  contactmap <- fread(file.path(dataset_dir, "processed_data", paste0("PDB_contactmap_", strsplit(basename(input_file), "\\.")[[1]][1], "_pseudomonomer_A", suffix, ".txt")))
  #Translate positions back to monomer positions
  contactmap[, Pos1 := (Pos1-1)%%nchar(aa_seq)+1]
  contactmap[, Pos2 := (Pos2-1)%%nchar(aa_seq)+1]
  #HAmin
  setkey(contactmap, HAmin)
  contactmap_HAmin <- contactmap[!duplicated(contactmap[,.(Pos1, Pos2)]), .(Pos1, Pos2, HAmin, HAmin_sd)]
  #scHAmin
  setkey(contactmap, scHAmin)
  contactmap_scHAmin <- contactmap[!duplicated(contactmap[,.(Pos1, Pos2)]), .(Pos1, Pos2, scHAmin, scHAmin_sd)]
  #CB
  setkey(contactmap, CB)
  contactmap_CB <- contactmap[!duplicated(contactmap[,.(Pos1, Pos2)]), .(Pos1, Pos2, CB, CB_sd)]
  #Merge
  setkey(contactmap, Pos1, Pos2)
  contactmap <- contactmap[!duplicated(contactmap[,.(Pos1, Pos2)]), .(Pos1, Pos2, WT_AA1, WT_AA2, chainids)]
  setkey(contactmap_HAmin, Pos1, Pos2)
  setkey(contactmap_scHAmin, Pos1, Pos2)
  setkey(contactmap_CB, Pos1, Pos2)
  contactmap <- contactmap[contactmap_HAmin,][contactmap_scHAmin,][contactmap_CB,][,.(Pos1, Pos2, WT_AA1, WT_AA2, chainids, HAmin, scHAmin, CB, HAmin_sd, scHAmin_sd, CB_sd)]
  #Save pairwise distance table
  pdb_filename = strsplit(strsplit(input_file,"/")[[1]][length(strsplit(input_file,"/")[[1]])],"\\.")[[1]][1]
  write.table(file = paste0(dataset_dir,'processed_data/PDB_contactmap_',pdb_filename,'_A',suffix,".txt",collapse = ""),
              x = contactmap,quote = F,row.names = F,col.names = T)
}

