#!/usr/bin/env Rscript
library(stringr)

#Create probability maps of pixels being in the nuclei
#Input are composite images in the well_location subdirectories
#Output are hf files in each well_location
#Call this script from eppec and from directory with well_location subdirectories
#that have a P_Reg subdirectories
#

#use command line arguments to identify the plate ID
# plate_path = commandArgs(trailingOnly=TRUE)
plate_path <- "/eppec/storage/groups/heiserlab/image_scratch/AU_I_L_008_01_1"
# dirs <- dir(path = plate_path, pattern = "[[:alnum:]]*_[[:digit:]]*",full.names = TRUE)[15]
dirs <- dir(path = plate_path, pattern = "C5_4|[D][12346]_[124]",full.names = TRUE)
#dirs <- dir(path = plate_path, pattern = "B1_4",full.names = TRUE)
res <- lapply(dirs, function(dir_path){
  message("classifying pixels in: ",dir_path)
  well_location <- str_remove(dir_path,".*/")
# Create pixel probability masks based on the high-contrast nuclear images and
# store them in the P_Reg directories
system(paste0('srun -c 8 --exclude=eppec-node6,eppec-node7 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --readonly=True --project=/graylab/share/dane/CellTracking/AU00801_Pixel_112918_plus_tracking.ilp --raw_data=',dir_path,'/Analysis/Composite_reg.tif --export_source="Probabilities" --output_format="hdf5" --output_filename_format=',dir_path,'/Analysis/Composite_reg_probs.h5'),
       wait=FALSE)
})

