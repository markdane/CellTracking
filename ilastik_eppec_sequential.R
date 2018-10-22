#!/usr/bin/env Rscript
library(stringr)
#Cell this script from a directory with Evos well_location subdirectories
#that have TxRed_Reg and GFP_Reg subdirectories
#####DEBUG skipping over first two subdirectories
pwd <- getwd()
dirs <- paste0(pwd,"/",dir(pattern = "[[:alnum:]]*_[[:digit:]]*"))[2:12]
#dir_path = "A1_1

res <- lapply(dirs, function(dir_path){
  message("processing images in: ",dir_path)
  well_location <- str_remove(dir_path,".*/")
  
  system(paste0('srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=/graylab/share/dane/CellTracking/HeLa_Pixel_Classifier.ilp --raw_data="',dir_path,'/TxRed_Reg/',well_location,'_TxRed_Scene1Interval*.tif" --output_format="hdf5" --output_filename_format=',dir_path,'/TxRed_Reg/Probabilities.h5'))
  
  system(paste0('srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=/graylab/share/dane/CellTracking/HeLa_Tracking.ilp --raw_data="',dir_path,'/TxRed_Reg/',well_location,'_TxRed_Scene1Interval*.tif" --prediction_maps=',dir_path,'/TxRed_Reg/Probabilities.h5 --export_source="Plugin" --export_plugin="CSV-Table"'))
  
  system(paste0('srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=/graylab/share/dane/CellTracking/HeLa_Tracking.ilp --raw_data="',dir_path,'/GFP_Reg/',well_location,'_GFP_Scene1Interval*.tif" --prediction_maps=',dir_path,'/TxRed_Reg/Probabilities.h5 --export_source="Plugin" --export_plugin="CSV-Table"'))
})

