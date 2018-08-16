##### extract secondary structure predictions from PSIPRED
SS_from_PSIPRED = function(input_file,
                           dataset_dir){
  
  secondary_structure = fread(input_file)[,.(Pos = V1,SS = V3)]
  write.table(paste0(dataset_dir,"processed_data/PSIPRED_secondary_structure.txt"),
              x = secondary_structure,quote = F,row.names = F,col.names = T)
  
}