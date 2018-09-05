
#################################################################
### convert pair distances into contact and distance matrices ###
#################################################################

contact_matrix_from_pairdistances <- function(input_file, 
  dataset_dir, 
  idx_start = 1,
  dist_type = c("scHAmin", "HAmin", "CB"),
  dist_cutoff = 4.5,
  plot = T){
  
  ### variables 
  # input_file: pairdistances data.table
  # dataset_dir: dataset directory, like "GB1/", 
  # idx_start: index of first position in contact matrix
  # dist_type: one of the following distance types: "scHAmin", "HAmin", "CB"
  # dist_cutoff: in Angstrom, used for determining the contact map
  # plot: whether to plot the heatmap (True, False)
  #### > will deposit contact maps (one for distance matrix and one for contact matrix) to dataset_dir/

  require(data.table)
  require(ggplot2)

  #Read pair distances file
  contactmap <- fread(input_file)
  #Build distance matrix
  all_positions <- unique(unlist(contactmap[,.(Pos1, Pos2)]))
  all_positions_names <- idx_start:(idx_start+length(all_positions)-1)
  num_positions <- length(all_positions)
  dist_mat <- matrix(nrow = num_positions, ncol = num_positions)
  for(i in all_positions){
    for(j in all_positions){
      dist_mat[i,j] <- as.numeric(contactmap[Pos1==i & Pos2==j,.SD,.SDcols = c(dist_type)])
    }
  }
  rownames(dist_mat) <- all_positions_names
  colnames(dist_mat) <- rownames(dist_mat)
  bin_mat<-dist_mat
  bin_mat[dist_mat<dist_cutoff]<-1
  bin_mat[dist_mat>=dist_cutoff]<-0
  bin_mat[dist_mat==0] <- NA
  #Plot heatmap matrices
  if(plot){
    tile_heatmap_wrapper(dist_mat[rev(1:dim(dist_mat)[1]),], file.path(dataset_dir, gsub(".txt$", "_distance_matrix.pdf", basename(input_file))), width=5, height=5, xlab = "Residue position", ylab = "Residue position", colour_clip=F, cluster='none', xaxis_size=10, yaxis_size=10, xaxis_angle=90, x_breaks = all_positions_names, y_breaks = all_positions_names)
    tile_heatmap_wrapper(bin_mat[rev(1:dim(bin_mat)[1]),], file.path(dataset_dir, gsub(".txt$", "_contact_matrix.pdf", basename(input_file))), width=5, height=5, xlab = "Residue position", ylab = "Residue position", colour_clip=F, cluster='none', xaxis_size=10, yaxis_size=10, xaxis_angle =90, x_breaks = all_positions_names, y_breaks = all_positions_names)
  }
  #Return
  return(list(distance_matrix = dist_mat, contact_matrix = bin_mat))
}

