#!/usr/bin/env Rscript
library(stringr)
#Call this script from a directory with well_location subdirectories
#that have R_Reg and G_Reg subdirectories
#####DEBUG skipping over first two subdirectories
pwd <- getwd()
dirs <- paste0(pwd,"/",dir(pattern = "[[:alnum:]]*_[[:digit:]]*"))[6:96]
#dir_path = "A3_1"

res <- lapply(dirs, function(dir_path){
  message("classifying pixels in nuclear images in: ",dir_path)
  well_location <- str_remove(dir_path,".*/")
# Create pixel probability masks based on the high-contrast nuclear images and
# store them in the R_Reg directories
  system(paste0('srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --readonly=True --project=/graylab/share/dane/CellTracking/HDF5_pixels.ilp --raw_data=',dir_path,'/P_Reg/Composite_reg.tif --output_format="hdf5" --output_filename_format=',dir_path,'/R_Reg/Probabilities.h5'),
         wait=FALSE)
})

###put a wait in here then the rest of it can run in parallel
#Use the nuclear probability masks in R_Reg to track the nuclei in the raw, registered images
#Export results to a csv table
# 
res <- lapply(dirs, function(dir_path){
  message("processing nuclear images in: ",dir_path)
  well_location <- str_remove(dir_path,".*/")
  system(paste0('srun -c 4 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --readonly=True --project=/graylab/share/dane/CellTracking/RGB_Tracking.ilp --raw_data="',dir_path,'/P_Reg/Composite_reg.tif" --prediction_maps=',dir_path,'/R_Reg/Probabilities.h5 --export_source="Plugin" --export_plugin="CSV-Table"'),
         wait=FALSE)
})

#   res <- lapply(dirs, function(dir_path){
#     message("processing reporter images in: ",dir_path)
#     well_location <- str_remove(dir_path,".*/")
#   system(paste0('srun -c 4 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --readonly=True --project=/graylab/share/dane/CellTracking/IncuCyte_tracking_classifier.ilp --raw_data="',dir_path,'/G_Reg/G_Unreg*.tif" --prediction_maps=',dir_path,'/R_Reg/Probabilities.h5 --export_source="Plugin" --export_plugin="CSV-Table"'),
#          wait=FALSE)
# })

