library(tidyverse)
library(EBImage)
library(rhdf5)
getGFPData <- function(redDirPath, GFPDirPath){
  red <- dir(redDirPath, pattern = "Scene1Interval[[:digit:]]*.tif", full.names = TRUE) %>%
    readImage()
  GFP <- dir(GFPDirPath, pattern =  "Scene1Interval[[:digit:]]*.tif", full.names = TRUE) %>%
    readImage()
  #Get nuclear probability 
  nuc_mask_dble <- dir(paste0(redDirPath,""), pattern = "Probabilities", full.names = TRUE) %>%
    h5read(name = "Beacon42") %>%
    .[2,,,]
  #convert probability masks to binary masks
  imageData(nuc_mask_dble)[imageData(nuc_mask_dble) > 0] <- 1
  # #invert binary masks
  # imageData(nuc_mask_dble) <- 1-imageData(nuc_mask_dble)
  
  #create cytoplasmic masks from nuc masks seeds
  ctmask = opening(GFP>0.1, makeBrush(5, shape='disc'))
  cell_mask = propagate(GFP, seeds=nuc_mask_dble, mask=ctmask)
  background_mask <- 1-cell_mask
  background <- GFP*background_mask
  #foo <- lapply(1:dim(background)[3], function(i){
  #  mean(imageData(background[,,i]))
 #})
  foo <- lapply(1:dim(background)[3], function(i){
    bg_pixel_int <- imageData(background[,,i])[imageData(background[,,i])>0]
    median(bg_pixel_int)
  })
}

# 
# write_csv(cyto_data, "A1_1/Cyto_data.csv")
