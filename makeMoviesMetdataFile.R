#!/usr/bin/env Rscript

#Label and tiem stamp movies from a set of directories
library(tidyverse)

#use command line arguments to identify the plate ID
#plateID = commandArgs(trailingOnly=TRUE)
data_dir <- "/eppec/storage/groups/heiserlab/image_scratch/"
plateID <- "LI_I_L_035_01_1"
plate_dir <- paste0(data_dir,plateID)
message("processing metadata for ",plateID)

if(!file.exists(paste0(plate_dir,"/Analysis/",plateID,"_movieMetadata.csv"))){
  #look in an Analysis directory for a csv file that holds the labels
  md <- dir(paste0(plate_dir,"/Analysis"),pattern = "csv", full.names = TRUE) %>%
    read_csv(col_types = cols(
      .default = col_character(),
      `Image-period` = col_character(),
      `CellLine-1-timeBegin` = col_character(),
      `CellLine-1-timeEnd` = col_character(),
      `Drug-1-conc` = col_character(),
      `Ligand-1-conc` = col_character(),
      `Ligand-1-timeBegin` = col_character(),
      `Ligand-2-conc` = col_character(),
      `Ligand-2-timeBegin` = col_character()
    ),
    na = " ")
  
  colnames(md) <- make.names(colnames(md))
  
  md <- md %>%
    select(Well, Image.period, Ligand.1, Ligand.2, Drug.1, Drug.1.conc, Drug.1.concUnit)
  md[is.na(md)] <- " "
  
  write_csv(md, paste0(plate_dir,"/Analysis/",plateID,"_movieMetadata.csv"))
  #Get the well_locations
}else {
  md <- read_csv(paste0(plate_dir,"/Analysis/",plateID,"_movieMetadata.csv"),
                 col_types = cols(
                   Well = col_character(),
                   Image.period = col_character(),
                   Ligand.1 = col_character(),
                   Ligand.2 = col_character(),
                   Drug.1 = col_character(),
                   Drug.1.conc = col_character(),
                   Drug.1.concUnit = col_character()
                 ),
                 na = " ")
}

#merge in the labels and time interval with the well_locations
dms <- dir(plate_dir, pattern = "[[:alnum:]]*_[[:digit:]]*") %>%
  tibble(Well_Location=.) %>%
  mutate(Well = str_remove(Well_Location, "_.*")) %>%
  left_join(md, by="Well")

#loop through all well_locations
foo <- lapply(seq_along(dms$Well_Location), function(i){
  subdirs <- dir(paste0(plate_dir,"/",dms$Well_Location[i]))
  lapply(subdirs[3], function(subdir){
   #run the imageJ macro to label images and wrote out AVI
    system(paste0("srun -c 1 ~/Fiji.app/ImageJ-linux64 --headless -macro /graylab/share/dane/CellTracking/label_movie.ijm ",
                  paste0('"',dms$Well_Location[i],"_",dms$Ligand.1[i],",",
                         paste0(plate_dir,"/",dms$Well_Location[i],","),
                         dms$Ligand.1[i],",",
                         dms$Ligand.2[i],",",
                         paste(dms$Drug.1[i],dms$Drug.1.conc[i],dms$Drug.1.concUnit[i],'"'))), wait = FALSE)
  })
})

