library(tidyverse)
library(EBImage)
library(rhdf5)

getGFPbackground <- function(redDirPath, GFPDirPath){
  red <- dir(redDirPath, pattern = "Scene1Interval[[:digit:]]*.tif", full.names = TRUE)%>%
    readImage()
  GFP_filenames <- dir(GFPDirPath, pattern =  "Scene1Interval[[:digit:]]*.tif", full.names = TRUE)
  GFP <- GFP_filenames %>%
    readImage()
  #Get nuclear probability 
  nuc_mask_dble <- dir(paste0(redDirPath,""), pattern = "Probabilities", full.names = TRUE) %>%
    h5read(name = "Beacon42") %>%
    .[2,,,]
  
  #convert nuclear probability masks to binary masks
  imageData(nuc_mask_dble)[imageData(nuc_mask_dble) > 0] <- 1

  #create cytoplasmic masks from nuc masks seeds
  ctmask = opening(GFP>0.1, makeBrush(5, shape='disc'))
  cell_mask = propagate(GFP, seeds=nuc_mask_dble, mask=ctmask)
  background_mask <- 1-cell_mask
  background <- GFP*background_mask

  bg_data <- lapply(1:dim(background)[3], function(i){
    bg_pixel_int <- imageData(background[,,i])[imageData(background[,,i])>0]
    bg_median <- median(bg_pixel_int)
    list(Filename = GFP_filenames[i], Background_median = bg_median)
  })
} %>%
  bind_rows()

#Calculate GFP  background and write out to files 
plate_dir <-  "/eppec/storage/groups/heiserlab/image_scratch/HE_E_L_049_01_1/"
well_locations <- dir(plate_dir, pattern = "[[:alnum:]]*_[[:digit:]]*", full.names = TRUE)
foo <- lapply(well_locations, function(well_location){
  res <- getGFPbackground(redDirPath = paste0(well_location,"/TxRed_Reg"),
                          GFPDirPath =  paste0(well_location,"/GFP_Reg"))
  write_csv(res, paste0(well_location,"/GFP_Reg/GFP_background.csv"))
})



