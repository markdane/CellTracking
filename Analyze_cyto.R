library(tidyverse)
library(EBImage)
getGFPData <- function(redDirPath, GFPDirPath){
  red <- dir(redDirPath, pattern = "Scene1Interval.*_TxRed.tif", full.names = TRUE) %>%
    readImage()
  GFP <- dir(GFPDirPath, pattern = "Scene1Interval.*_GFP.tif", full.names = TRUE) %>%
    readImage()
  nuc_mask_dble <- dir(paste0(redDirPath,""), pattern = ":Users", full.names = TRUE) %>%
    readImage()
  imageData(nuc_mask_dble) <- imageData(nuc_mask_dble)*2^16

  nuc_mask_labels <- lapply(1:dim(nuc_mask_dble)[3], function(i){
    apply(imageData(nuc_mask_dble[,,i]), 2, as.integer)
  })
  
  nuc_mask <- Image(unlist(nuc_mask_labels), dim = c(dim(nuc_mask_labels[[1]]),length(nuc_mask_labels)))
  # cells_Hela <- rgbImage(green = 1.5*cel_Hela, blue = nuc_Hela)
  # display(cells_Hela, all = TRUE)
  # display(colorLabels(nmask_Hela), all=TRUE)
  
  ctmask = opening(GFP>0.1, makeBrush(5, shape='disc'))
  cell_mask = propagate(GFP, seeds=nuc_mask, mask=ctmask)
  cyto_mask <- cell_mask - nuc_mask
  
  #display(colorLabels(red), all=TRUE)
  #display(colorLabels(nuc_mask_dble), all=TRUE)
  #display(colorLabels(nuc_mask), all=TRUE)
  #display(colorLabels(cyto_mask), all=TRUE)
  #display(colorLabels(cell_mask), all=TRUE)
  # segmented = paintObjects(cmask, cells_Hela, col='#ff00ff')
  # segmented = paintObjects(nmask_Hela, segmented, col='#ffff00')
  # display(segmented, all=TRUE)
  GFP_data <- lapply(1:(dim(cell_mask)[3]), function(frame){
    cell_data <- computeFeatures.basic(cell_mask[,,frame], GFP[,,frame]) %>%
      as.tibble() %>%
      rename_at(names(.), funs( paste0("Cell_",.)))
    cell_data <- tibble(frame = frame-1 , trackId_nuclear = 1:nrow(cell_data)) %>%
      bind_cols(cell_data)
    nuclear_data <- computeFeatures.basic(nuc_mask[,,frame], GFP[,,frame]) %>%
      as.tibble() %>%
      rename_at(names(.), funs( paste0("Nuclei_",.) ) ) %>%
      mutate(frame = frame -1,
             trackId_nuclear = 1:nrow(.))
    cyto_data <- computeFeatures.basic(cyto_mask[,,frame], GFP[,,frame]) %>%
      as.tibble() %>%
      rename_at(names(.), funs( paste0("Cyto_",.) ) )%>%
      mutate(frame = frame -1,
             trackId_nuclear = 1:nrow(.))
    all_data <- cell_data %>%
      left_join(nuclear_data, by = c("frame", "trackId_nuclear")) %>%
      left_join(cyto_data, by = c("frame", "trackId_nuclear"))
  }) %>%
    bind_rows()
}

cyto_data <- getGFPData(redDirPath = "Beacon-42/Red_Reg",GFPDirPath = "Beacon-42/GFP_Reg")
write_csv(cyto_data, "Beacon-42/Cyto_data.csv")
