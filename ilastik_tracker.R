#!/usr/bin/env Rscript
library(stringr)

#Create csv tables of feature data for tracked cells
#Input are probability masks in the well_location subdirectories
#Output are csv files in each well_location
#Call this script from eppec with a plate path command line argument

# plate_path = commandArgs(trailingOnly=TRUE)
plate_path <- "/eppec/storage/groups/heiserlab/image_scratch/AU_I_L_008_01_1"
#dirs <- dir(path = plate_path, pattern = "[[:alnum:]]*_[[:digit:]]*",full.names = TRUE)[c(17:19)]
dirs <- dir(path = plate_path, pattern = "C4_2|D4_4",full.names = TRUE)
res <- lapply(dirs, function(dir_path){
  message("tracking cells in: ",dir_path)
# Create csv files baed on nuclei probability masks
  system(paste0('srun -c 8  --exclude=eppec-node4 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --readonly=True --project=/graylab/share/dane/CellTracking/Tracking_AU00801_v4.ilp --raw_data=',dir_path,'/Analysis/Composite_reg.tif --prediction_maps=',dir_path,'/Analysis/Composite_reg_probs.h5 --export_source="Plugin" --export_plugin="CSV-Table"'),
         wait=FALSE)
})

