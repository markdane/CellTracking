#!/usr/bin/env Rscript
library(stringr)

#Create csv table of feature data for tracked cells
#Input are probability masks in the well_location subdirectories
#Output are csv files in each well_location
#Call this script from eppec with a plate path command line argument

#use command line arguments to identify the plate ID
#plate_path = commandArgs(trailingOnly=TRUE)
plate_path <- "/eppec/storage/groups/heiserlab/image_scratch/LI_I_L_035_01_1"
dirs <- dir(path = plate_path, pattern = "[[:alnum:]]*_[[:digit:]]*",full.names = TRUE)

res <- lapply(dirs, function(dir_path){
  message("tracking cells in: ",dir_path)
  well_location <- str_remove(dir_path,".*/")
# Create pixel probability masks based on the high-contrast nuclear images and
# store them in the P_Reg directories
  system(paste0('srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --readonly=True --project=/graylab/share/dane/CellTracking/RGB_Tracking.ilp --raw_data="',dir_path,'/Composite_reg.tif" --prediction_maps=',dir_path,'/Probabilities.h5 --export_source="Plugin" --export_plugin="CSV-Table"'),
         wait=FALSE)
})

