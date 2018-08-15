#this is the script to initialize R on the cluster and call the actual XPLOR modeling function

require(data.table)
require(optparse)

#read in varlist file deposited in the copied folder structure
option_list = list(
  make_option(opt_str = c("-v", "--varlist"), type="character", default=NULL,
              help="varlist pointer", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

file_pointer = opt$varlist
load(file = file_pointer)

###### source all XPLOR modeling functions ###### 
source(paste0(varlist$cluster_dir,varlist$protein,"/XPLOR_modeling_functions_v2.R"))

###### start XPLOR modeling
XPLOR_modeling(varlist)