#Structure files and directories for processing
library(tidyverse)
raw_data_path <-  "/graylab/share/grossse/AU00701/"
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

createDirs <- function(dir_names, raw_data_path){
  foo <- lapply(dir_names, function(dir_name){
    full_dir_name <- paste0(raw_data_path,dir_name)
    if(!dir.exists(full_dir_name)) dir.create(full_dir_name)
    return(full_dir_name)
  })
}

copyFiles <- function(well_location, filename){
  to_filename <- str_extract(filename, ".*/" ) %>%
    paste0(well_location, "/",str_remove(filename, ".*/" ))
  #convert this to system(paste("mv",filename, to_filename))
  file.copy(filename, to_filename)
}

#create directories for each well + location
foo <- files %>%
  select(Well_Location) %>%
  distinct() %>%
  unlist() %>%
  createDirs(raw_data_path)

#Move files to the correct directory
foo <- files %>%
  select(Well_Location, Full_filename) %>%
  group_by(Well_Location) %>%
  mutate(foo = copyFiles(Well_Location, Full_filename))

