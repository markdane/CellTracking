#!/usr/bin/env Rscript

#Label and tiem stamp movies from a set of directories
library(tidyverse)
library(magick)
library(EBImage)

#use command line arguments to identify the plate ID
plateID = commandArgs(trailingOnly=TRUE)
plateID <- "/eppec/storage/groups/heiserlab/image_scratch/LI_I_L_035_01_1/"
message("processing images in ",plateID)

#look in an Analysis directory for a csv file that holds the labels
md <- dir(paste0(plateID,"Analysis"),pattern = "csv", full.names = TRUE) %>%
  read_csv()

colnames(md) <- make.names(colnames(md))

md <- md %>%
  select(Well, Image.period, Ligand.1, Ligand.2, Drug.1, Drug.1.conc, Drug.1.concUnit)

#Get the well_locations
#merge in the labels and time interval with the well_locations

dms <- dir(plateID, pattern = "[[:alnum:]]*_[[:digit:]]*") %>%
  tibble(Well_Location=.) %>%
  mutate(Well = str_remove(Well_Location, "_.*")) %>%
  left_join(md, by="Well")


#loop through all well_locations
foo <- lapply(dms$Well_Location[1], function(dm){
  subdirs <- dir(paste0(plateID,dm))
  lapply(subdirs[3], function(subdir){
    #read in the images
    img_files <- dir(paste0(plateID,dm,"/",subdir), pattern = "tif", full.names = TRUE)[1:2] 
    foo <- lapply(img_files, function(imgFn){
      img <- readImage(imgFn,all = FALSE)
      display(img, method="raster")
      text(x = 1200, y=1000, labels = "OSM", adj = c(0,1), col = "white", cex = 2)
      filename = "test"
      dev.print(jpeg, filename = filename , width = dim(img)[1], height = dim(img)[2])
    })

      read()
    #add the labels
    
    foo <- image_annotate(images,text = "skdjfsk")
    #add the timestamps
    #write the movies to the Analysis directory
  })


})

