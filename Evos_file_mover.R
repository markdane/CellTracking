#Structure files and directories for processing
library(tidyverse)

#Execute script from the plate level
#Will create and copy png or tiff files from Beacon structure to well location structure
#Requires meatdata file with Well and Beacon columns
raw_data_path <-  getwd()
plateDir <- str_remove(raw_data_path, ".*/")

#Read in the metadata file
metadata <- dir("Analysis/", pattern = "metadata", full.names = TRUE) %>%
  read_csv() %>%
  select(Well, Beacon) %>%
  arrange(Beacon) %>%
  group_by(Well) %>%
  mutate(Location = row_number(),
         Well_Location = paste0(Well,"_",Location))
if(nrow(metadata)==0) stop("No metadata file found in ", getwd(),"/Analysis/" )

#Get a tibble of all files in the dataset and split to create metadata
files<- tibble(Full_filename=dir(pattern = "GFP|TxRed", full.names = TRUE, recursive = TRUE)) %>%
  mutate(Filename = str_remove(Full_filename, ".*/"),
         Plate_Dir = plateDir,
         Channel = str_remove(Filename,"[[:alnum:]]*_{1}"),
         Channel = str_remove(Channel,".png|.tif"),
         Beacon = str_extract(Full_filename,"Beacon-[[:digit:]]*"),
         Beacon = str_remove(Beacon, "Beacon-"),
         Beacon = as.integer(Beacon),
         Index = str_remove(Filename, "_.*"),
         Filetype = str_extract(Filename,".png|.tif|.tiff")) %>%
  left_join(metadata,by = "Beacon") %>%
  mutate(
    To_Filename = paste0(Well_Location, "_", Channel, "_",Index, Filetype))

createDirs <- function(dir_names, raw_data_path){
  foo <- lapply(dir_names, function(dir_name){
    full_dir_name <- paste0(raw_data_path,"/",dir_name)
    if(!dir.exists(full_dir_name)) dir.create(full_dir_name)
    return(full_dir_name)
  })
}

copyFiles <- function(filename, well_location, to_filename){
  to_filepath <- paste0(well_location, "/", to_filename)
  #convert this to system(paste("mv",filename, to_filename))
  file.copy(filename, to_filepath)
}

#create directories for each well + location
foo <- files %>%
  select(Well_Location) %>%
  distinct() %>%
  unlist() %>%
  createDirs(raw_data_path)

#Move files to the correct directory
  foo <- files %>%
    select(Well_Location, Full_filename, To_Filename) %>%
    group_by(Well_Location) %>%
  mutate(foo = copyFiles(Full_filename, Well_Location, To_Filename))

