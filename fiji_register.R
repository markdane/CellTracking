#!/usr/bin/env Rscript
library(stringr)

#register images
#Input are raw image stacks in R_Unreg, G_Unreg and P_Unreg well_location subdirectories
#Ouptut are composite images in the well_location subdirectories
#Call this script from eppec with a plate path command line argument

#use command line arguments to identify the plate ID
plate_path = commandArgs(trailingOnly=TRUE)
#plate_path <- "/eppec/storage/groups/heiserlab/image_scratch/LI_I_L_035_01_1"
dirs <- dir(path = plate_path, pattern = "[[:alnum:]]*_[[:digit:]]*",full.names = TRUE)

res <- lapply(dirs, function(dir_path){
  message("registering images in ",dir_path)
# Call a fiji macro on each well_location directory
# store composite image results in the P_Reg directories
  system(paste0('srun -c 2 ~/Fiji.app/ImageJ-linux64 --headless -macro /graylab/share/dane/CellTracking/RGB_IC_Arg.ijm ',dir_path),
         wait=TRUE)
})

