#!/usr/bin/env Rscript
library(stringr)

#Create probability maps of pixels being in the nuclei
#Input are composite images in the well_location subdirectories
#Output are hf files in each well_location
#Call this script from eppec and from directory with well_location subdirectories
#that have a P_Reg subdirectories
#

#use command line arguments to identify the plate ID
plate_path = commandArgs(trailingOnly=TRUE)
plate_path <- "/eppec/storage/groups/heiserlab/image_scratch/LI_I_L_035_01_1"
dirs <- dir(path = plate_path, pattern = "[[:alnum:]]*_[[:digit:]]*",full.names = TRUE)

res <- lapply(dirs, function(dir_path){
  message("classifying pixels in: ",dir_path)
  well_location <- str_remove(dir_path,".*/")
# Create pixel probability masks based on the high-contrast nuclear images and
# store them in the P_Reg directories
system(paste0('srun -c 4 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --readonly=True --project=/graylab/share/dane/CellTracking/HDF5_pixels.ilp --raw_data=',dir_path,'/Analysis/Composite_reg.tif --output_format="hdf5" --output_filename_format=',dir_path,'/Analysis/Probabilities.h5'),
       wait=FALSE)
})

