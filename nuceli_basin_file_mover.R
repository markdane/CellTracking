#Structure files and directories for processing
library(tidyverse)
raw_data_path <-  "/graylab/share/dataxchange/LincsWorkflow/"
tiff_data_path <- "/graylab/share/dane/MEP-LINCS/DevelopmentExperiments/MEMA_Analysis/ImageCache/HCC1143_COL1/"
#Start in top level subdrectory
#Get a tibble of the files in the dataset and split to create metadata
files<- tibble(Full_filename=dir(paste0(raw_data_path,"HCC1143_COL1 - Segmentation") ,pattern = "Nuclei Basins", full.names = TRUE)) %>%
  mutate(Filename = str_remove(Full_filename, ".*/"),
         ImageID = str_remove(Filename, " .*"))

basin_dirs <- dir(raw_data_path, pattern = "Group._Basins")
foo <- lapply(basin_dirs, function(basin_dir){
  group <- str_extract(basin_dir, "Group[[:digit:]]") %>%
    str_remove("Group") %>%
    as.integer()
  lapply(((group-1)*1000+1):((group*1000)), function(i){
    file.copy(files$Full_filename[i], paste0(raw_data_path,basin_dir,"/",files$Filename[i]))
    #foo <-  c(paste0(tiff_data_path,files$ImageID[i],.tif), paste0(raw_data_path,basin_dir,"/",files$Filename[i]))
  })
})

tiff_dirs <- dir(raw_data_path, pattern = "Group.$")
foo <- lapply(tiff_dirs, function(tiff_dir){
  group <- str_extract(tiff_dir, "Group[[:digit:]]") %>%
    str_remove("Group") %>%
    as.integer()
  lapply(((group-1)*1000+1):((group*1000)), function(i){
    file.copy(paste0(tiff_data_path,files$ImageID[i],".tif"), paste0(raw_data_path,tiff_dir,"/",files$ImageID[i],".tif"))
    #foo <-  c(paste0(tiff_data_path,files$ImageID[i],".tif"), paste0(raw_data_path,tiff_dir,"/",files$ImageID[i],".tif"))
  })
})

