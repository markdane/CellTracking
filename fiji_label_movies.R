#!/usr/bin/env Rscript

#Label movies in plate subdirectories
library(tidyverse)

#use command line arguments to identify the plate ID
plate_path = commandArgs(trailingOnly=TRUE)
#plate_path <- "/eppec/storage/groups/heiserlab/image_scratch/LI_I_L_035_01_1"
#plate_path <- "/Users/dane/Documents/image_scratch/LI_I_L_035_01_1"
plateID <- str_remove(plate_path,".*/")
message("processing metadata for ",plateID)

#look in an Analysis directory for a csv file that holds the labels
  md <- dir(paste0(plate_path,"/Analysis"),pattern = "csv", full.names = TRUE) %>%
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
  
#merge in the labels and time interval with the well_locations
dms <- dir(plate_path, pattern = "[[:alnum:]]*_[[:digit:]]*") %>%
  tibble(Well_Location=.) %>%
  mutate(Well = str_remove(Well_Location, "_.*")) %>%
  left_join(md, by="Well")

#loop through all well_locations
foo <- lapply(seq_along(dms$Well_Location), function(i){
   #run the imageJ macro to label images and wrote out AVI
    system(paste0("srun -c 1 ~/Fiji.app/ImageJ-linux64 --headless -macro /graylab/share/dane/CellTracking/label_movie.ijm ",
                  paste0('"',dms$Well_Location[i],"_",dms$Ligand.1[i],",",
                         paste0(plate_path,"/",dms$Well_Location[i],","),
                         dms$Ligand.1[i],",",
                         dms$Ligand.2[i],",",
                         paste(dms$Drug.1[i],dms$Drug.1.conc[i],dms$Drug.1.concUnit[i],'"'))), wait = FALSE)
})

