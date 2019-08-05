#Structure files and directories for processing
library(tidyverse)
raw_data_path <-  "/data/share/grossse/AU_I_L_007_01_1/AU_I_L_007_01_1_IMAGES"
new_raw_data_path <- "/eppec/storage/groups/heiserlab/image_scratch/AU_I_L_007_01_1"
#Start in top level subdrectory
#Get a tibble of the files in the dataset and split to create metadata
files<- tibble(Full_filename=dir(raw_data_path,pattern = "tif", full.names = TRUE)) %>%
  mutate(Filename = str_remove(Full_filename, ".*/"),
         Barcode = str_remove(Filename, "_.*"),
         Channel = str_remove(Filename,"[[:alnum:]]*_{1}"),
         Channel = str_remove(Channel,"_.*"),
         Well = str_remove(Filename,"[[:alnum:]]*_[[:alnum:]]*_"),
         Well = str_remove(Well, "_.*"),
         Location = str_remove(Filename,"[[:alnum:]]*_[[:alnum:]]*_[[:alnum:]]*_"),
         Location = str_remove(Location, "_.*"),
         Location = as.integer(Location),
         Well_Location = paste(Well, Location, sep="_"),
         Slice = str_remove(Filename, ".*_"),
         Slice = str_remove(Slice, ".tif"))

#beginning of a cleaner way to decode the metadata
#foo <- lapply(beacon_list, str_split, pattern="_")

createDirs <- function(dir_names){
  #create well _location directories first
  well_location_dirs <- dir_names %>%
    select(New_dir_name) %>%
    distinct() %>%
    unlist()
  dir_success <- lapply(well_location_dirs, function(well_location_dir){
    if(!dir.exists(well_location_dir)) dir.create(well_location_dir)
  })
  
  unreg_dirs <- dir_names %>%
    mutate(New_channel_dir_name = paste0(New_dir_name,"/",Channel,"_Unreg")) %>%
    select(New_channel_dir_name) %>%
    unlist()
  
  dir_success <- lapply(unreg_dirs, function(unreg_dir){
    if(!dir.exists(unreg_dir)) dir.create(unreg_dir)
  })
  
}

copyFiles <- function(filename, new_file_path){
  #convert this to system(paste("mv",filename, to_filename))
  foo <- lapply(seq_along(filename), function(i){
    file.copy(filename[i], new_file_path[i])
  })
}

#create directories for each well + location
foo <- files %>%
  select(Well_Location, Channel) %>%
  distinct() %>%
  mutate(New_dir_name = paste0(new_raw_data_path,"/",Well_Location)) %>%
  createDirs()

#Move files to the correct directory
foo <- files %>%
  select(Well_Location, Channel, Full_filename) %>%
  mutate(New_file_path = paste0(new_raw_data_path,"/",Well_Location,"/",Channel,"_Unreg")) %>%
  mutate(foo = copyFiles(Full_filename, New_file_path))
