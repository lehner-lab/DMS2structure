create_directory_structure = function(directory) {
  ### create directory and subdirectories for a dataset
  #e.g. directory = "GB1/"
  system(command = paste0("mkdir -p ",directory), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"dataset/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"dataset/PDB/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"dataset/PSIPRED/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"processed_data/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"results/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"results/epistasis/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"results/preprocessing/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"results/secondary_structure/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"results/secondary_structure/temp_plots/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"results/tertiary_contacts/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"results/deepcontact/")), wait = T)
  system(command = paste0("mkdir -p ",paste0(directory,"results/XPLOR/")), wait = T)

}