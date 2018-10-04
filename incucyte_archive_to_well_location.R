#Structure files and directories for processing
library(tidyverse)
raw_data_path <-  "/graylab/share/grossse/AU00601_INCUCYTE"
#Start in top level subdrectory
#Get a tibble of the files in the dataset and split to create metadata
files <- str_split(dir(raw_data_path,pattern = "tif", full.names = TRUE,recursive = TRUE), "/", simplify = TRUE) %>%
  as.tibble() %>%
  rename(Year_Month = V8,
         Day = V9,
         Hour_Minute = V10,
         VesselID = V11,
         Filename = V12) %>%
  mutate(Dir_path = paste(V1,V2,V3,V4,V5,V6,V7, sep="/"),
         Well = str_remove(Filename,"-.*"),
         Location = str_extract(Filename, "-.*-"),
         Location = str_extract(Filename,"[[:digit:]*]"),
         New_filename = str_remove(Filename,".tif"),
         New_filename = paste0(New_filename,"_",Year_Month,"_",Day,"_",Hour_Minute,".tif"),
         New_filename = str_replace_all(New_filename,"-","_"),
         Well_Location = paste(Well, Location, sep="_"),
         Full_filename = dir(raw_data_path,pattern = "tif", full.names = TRUE,recursive = TRUE)) %>%
  select(-V1, -V2, -V3, -V4, -V5, -V6, -V7)

createDirs <- function(dir_names, raw_data_path){
  foo <- lapply(dir_names, function(dir_name){
    full_dir_name <- paste0(raw_data_path,"/",dir_name)
    if(!dir.exists(full_dir_name)) dir.create(full_dir_name)
    return(full_dir_name)
  })
}

copyFiles <- function(well_location, new_filename, filename){
  to_filename <- paste0(raw_data_path,"/",well_location,"/",new_filename)
  #convert this to system(paste("mv",filename, to_filename))
  system(paste("mv",filename, to_filename))
}

#create directories for each well + location
foo <- files[1,] %>%
  select(Well_Location) %>%
  distinct() %>%
  unlist() %>%
  createDirs(raw_data_path)

#Move files to the correct directory
foo <- files[1,] %>%
  select(Well_Location, New_filename, Full_filename) %>%
  group_by(Well_Location) %>%
  mutate(foo = copyFiles(Well_Location, New_filename, Full_filename))

