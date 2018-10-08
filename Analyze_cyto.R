library(tidyverse)
library(EBImage)
getGFPData <- function(redDirPath, GFPDirPath){
  red <- dir(redDirPath, pattern = "Scene1Interval.*_TxRed.tif", full.names = TRUE) %>%
    readImage()
  GFP <- dir(GFPDirPath, pattern = "Scene1Interval.*_GFP.tif", full.names = TRUE) %>%
    readImage()
  nuc_mask <- dir(redDirPath, pattern = "Scene1Interval0________Tracking-Result_.*", full.names = TRUE) %>%
    readImage() %>%
    bwlabel()
  
  # cells_Hela <- rgbImage(green = 1.5*cel_Hela, blue = nuc_Hela)
  # display(cells_Hela, all = TRUE)
  # display(colorLabels(nmask_Hela), all=TRUE)
  
  ctmask = opening(GFP>0.1, makeBrush(5, shape='disc'))
  cell_mask = propagate(GFP, seeds=nuc_mask, mask=ctmask)
  cyto_mask <- cell_mask - nuc_mask
  
  #display(colorLabels(cyto_mask), all=TRUE)
  # segmented = paintObjects(cmask, cells_Hela, col='#ff00ff')
  # segmented = paintObjects(nmask_Hela, segmented, col='#ffff00')
  # display(segmented, all=TRUE)
  GFP_data <- lapply(1:dim(cell_mask)[3], function(slice){
    cell_data <-     computeFeatures.basic(cell_mask[,,slice], GFP[,,slice]) %>%
      as.tibble() %>%
      rename_at(names(.), funs( paste0("Cell_",.)))
    cell_data <- tibble(Slice = slice, TrackId = 1:nrow(cell_data)) %>%
      bind_cols(cell_data)
    nuclear_data <- computeFeatures.basic(nuc_mask[,,slice], GFP[,,slice]) %>%
      as.tibble() %>%
      rename_at(names(.), funs( paste0("Nuclei_",.) ) ) %>%
      mutate(Slice = slice,
             TrackId = 1:nrow(.))
    cyto_data <- computeFeatures.basic(cyto_mask[,,slice], GFP[,,slice]) %>%
      as.tibble() %>%
      rename_at(names(.), funs( paste0("Cyto_",.) ) )%>%
      mutate(Slice = slice,
             TrackId = 1:nrow(.))
    all_data <- cell_data %>%
      left_join(nuclear_data, by = c("Slice", "TrackId")) %>%
      left_join(cyto_data, by = c("Slice", "TrackId"))
  }) %>%
    bind_rows()
}

all_data <- getGFPData(redDirPath = "Beacon-42/Red_Reg",GFPDirPath = "Beacon-42/GFP_Reg")
