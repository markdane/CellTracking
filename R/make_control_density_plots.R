library(tidyverse)
library(viridis)

read_plate_l1 <- function(plateID){
  l1_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_level_1.csv")
    print(paste("reading level 1 control data for ",plateID))
    l1_control_data <- map_dfr(dir(paste0(data_path,plateID,"/Analysis/",pipeline_name),pattern = paste0(plateID,".*(A1|A3|A5)_.*_level_1.csv"), full.names = TRUE), read_csv, show_col_types = FALSE, lazy = FALSE,
                       col_types =cols(.default = "d",
                                       well = col_character(),
                                       time_slice = col_character(),
                                       is_parent= col_logical(),
                                       plateID = col_character(),
                                       Well = col_character(),
                                       Cellline = col_character(),
                                       Ligand1 = col_character(),
                                       Ligand1ConcentrationUnits = col_character(),
                                       Ligand2 = col_character(),
                                       Ligand2ConcentrationUnits = col_character(),
                                       Ligand3= col_character(),
                                       Ligand3ConcentrationUnits = col_character(),
                                       ECMP1 = col_character(),     
                                       ECMP1ConcentrationUnits = col_character(),
                                       Drug1 = col_character(),
                                       Drug1Concentration = col_character(),
                                       Drug1ConcentrationUnits = col_character(),
                                       Drug2 = col_character(),
                                       Drug2Concentration = col_character(),
                                       Drug2ConcentrationUnits = col_character(),
                                       row = col_character(),
                                       Endpoint_red = col_character(),
                                       Endpoint_green = col_character(),
                                       treatment = col_character()))  %>%
      mutate( treatment =  paste(Drug1, Drug1Concentration,Drug2, Drug2Concentration, sep = "_"),
              treatment = str_remove_all(treatment, "_none|_0"))

  return(list(l1_control_data = l1_control_data))
}
  
create_density_plots <- function(plateID){
  #browser()
  l1_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_level_1.csv")
  
  df <- datasets[[plateID]][["l1_control_data"]] 
  
  df_select <- df %>%
    mutate(hours = elapsed_minutes/60) %>%
    select(Cell_CC_mean_intensity_ratio, treatment, hours) %>%
    filter(hours %in% c(0, 24, 48, 72, 96)) %>%
    drop_na()
  
  cols <- c("control"="royalblue", "vehicle"="royalblue", "Untreated"="royalblue")
  
  plateID <- unique(df$plateID)
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_",pipeline_name,"_contol_density_plots.pdf"))
  
  p_untreated_densities <- ggplot(df_select, aes(Cell_CC_mean_intensity_ratio, fill = treatment, color = treatment)) +
    geom_density(alpha = .2)+
    #geom_vline(xintercept = .8)+
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    facet_wrap(~hours) +
    guides(fill = "none", color = "none") +
    theme_bw()
  
  print(p_untreated_densities)
  
  res <- dev.off()
}


###########

data_path <-  "/home/exacloud/gscratch/HeiserLab/images/"
pipeline_name <- "CtcK"
plateIDs <- c("HC01301" = "HC01301","HC01401" = "HC01401")


datasets <- map(plateIDs, read_plate_l1)
res <- map(plateIDs, create_density_plots)
