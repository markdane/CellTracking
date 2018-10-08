#Structure files and directories for processing
library(tidyverse)
#raw_data_path <-  "/graylab/share/grossse/AU_I_L_006_01_1/AU_I_L_006_01_1/"
#Start in top level subdrectory
#Get a tibble of the files in the dataset and split to create metadata
files <- str_split(dir(pattern = "tif", full.names = TRUE,recursive = TRUE), "/", simplify = TRUE) %>%
  as.tibble() %>%
  rename(Year_Month = V4,
         Day = V5,
         Hour_Minute = V6,
         VesselID = V7,
         Filename = V8) %>%
  mutate(Dir_path = paste(V1,V2,V3, sep="/"),
         Well = str_remove(Filename,"-.*"),
         Location = str_extract(Filename, "-.*-"),
         Location = str_extract(Location,"[[:digit:]*]"),
         New_filename = str_remove(Filename,".tif"),
         New_filename = paste0(New_filename,"_",Year_Month,"_",Day,"_",Hour_Minute,".tif"),
         New_filename = str_replace_all(New_filename,"-","_"),
         Well_Location = paste(Well, Location, sep="_"),
         Full_filename = dir(pattern = "tif", full.names = TRUE,recursive = TRUE)) %>%
  select(-V1, -V2, -V3)

createDirs <- function(dir_names){
  foo <- lapply(dir_names, function(dir_name){
    full_dir_name <- paste0(dir_name)
    if(!dir.exists(full_dir_name)) dir.create(full_dir_name)
    return(full_dir_name)
  })
}

copyFiles <- function(well_location, new_filename, filename){
  to_filename <- paste0(well_location,"/",new_filename)
  #convert this to system(paste("mv",filename, to_filename))
  foo <- lapply(1:length(to_filename), function(i){
    system(paste("mv",filename[i], to_filename[i]))
    #system(paste("cp",filename[i], to_filename[i]))
    
  })
}

#create directories for each well + location
foo <- files %>%
  select(Well_Location) %>%
  distinct() %>%
  unlist() %>%
  createDirs()

#Move files to the correct directory
foo <- files %>%
  select(Well_Location, New_filename, Full_filename) %>%
  group_by(Well_Location) %>%
  mutate(foo = copyFiles(Well_Location, New_filename, Full_filename))

