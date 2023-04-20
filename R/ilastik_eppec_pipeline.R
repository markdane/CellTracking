#!/usr/bin/env Rscript

#Cell this script from a directory with Evos well_location subdirectories
#that have TxRed_Reg and GFP_Reg subdirectories
library(optparse)

#Get the path to the plateID from the command line

# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  parser <- OptionParser(usage = "%prog [options] WELL_LOCATION_DIR")
  arguments <- parse_args(parser, positional_arguments = 1)
}

#Specify the command line options
cl <- getCommandLineArgs()
dir_path <- cl$args[1]
#dir_path = "A1_1
message("processing images in: ",dir_path)

system(paste0('srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=/graylab/share/dane/CellTracking/HeLa_Pixel_Classifier.ilp --raw_data="./',dir_path,'/TxRed_Reg/',dir_path,'_TxRed_Scene1Interval*.tif" --output_format="hdf5" --output_filename_format=./',dir_path,'/TxRed_Reg/Probabilities.h5'))
 
system(paste0('srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=/graylab/share/dane/CellTracking/HeLa_Tracking.ilp --raw_data="./',dir_path,'/TxRed_Reg/',dir_path,'_TxRed_Scene1Interval*.tif" --prediction_maps=./',dir_path,'/TxRed_Reg/Probabilities.h5  --export_source="Plugin" --export_plugin="CSV-Table"'))
  
  system(paste0('srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=/graylab/share/dane/CellTracking/HeLa_Tracking.ilp --raw_data="./',dir_path,'/GFP_Reg/',dir_path,'_GFP_Scene1Interval*.tif" --prediction_maps=./',dir_path,'/TxRed_Reg/Probabilities.h5  --export_source="Plugin" --export_plugin="CSV-Table"'))
