library(tidyverse)

calc_distance <- function(x, y){
  if( length(x) == 1) return(d=0)
  d <- numeric(length(x)) #intialize a vector with 0s 
  for(i in 2:length(x)) { #start distances in 2nd element
    d[i] <- sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2)
  }
  return(round(d, digits = 2))
}

calc_direction <- function(x, y){
  if( length(x) == 1) return(d=0)
  d <- numeric(length(x)) #intialize a vector with 0s 
  for(i in 2:length(x)) { #start directions in 2nd element
    d[i] <- atan2(y[i]-y[i-1], x[i]-x[i-1])*180/pi
  }
  return(round(d, digits = 2))
}

#recursive function to return lineage
process_tracks_file <- function(tracks_filename){
  
  create_lineage <- function(x){
    mother <- subset(tracks, label==x)$parent
    if(mother==0) return(c(x,0))
    #browser()
    #print(x, mother)
    c(x,create_lineage(mother))
  }
  
  tracks <- read_csv(tracks_filename,
                     col_types = cols(
                       label = col_double(),
                       begins = col_double(),
                       ends = col_double(),
                       parent = col_double(),
                       length = col_double(),
                       is_parent = col_logical()
                     ))
  
  #Create a lineage column
  #get the lineage tree and name it by the second to last cell
  
  tracks_lineages <- tracks %>%
    mutate(lineage_labels = map(label, create_lineage),
           lineage = lapply(lineage_labels, dplyr::nth,-2),
           lineage = unlist(lineage)) %>%
    select(plateID, well, field,  label, lineage)
  
  #save this version of the filtered tracks file with lineage info
  res <- write_csv(tracks_lineages, str_replace(tracks_filename, "tracks", "tracks_lineages"))
  return(tracks_lineages)
}

read_plate_l1 <- function(plateID){
  l1_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_level_1.csv")
  
  print(paste("creating new level 1 file for ",plateID))
  l1_data <- map_dfr(dir(paste0(data_path,plateID,"/Analysis/",pipeline_name),pattern = paste0(plateID,"_.*_level_1.csv"), full.names = TRUE), read_csv, show_col_types = FALSE, lazy = FALSE,
                     col_types =cols(.default = "d",
                                     Nuclei_CKn_CC_mean_intensity = col_double(),
                                     well = col_character(),
                                     time_slice = col_character(),
                                     cell_cycle_state = col_character(),
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
                                     'incubation _media' = col_character(),
                                     treatment = col_character())) %>%
    mutate(plateID = plateID,
           day = str_remove(time_slice, "d.*"),
           day = as.integer(day),
           hour = str_extract(time_slice, "[0-9]*h"),
           hour = str_remove(hour, "h"),
           hour = as.integer(hour),
           minute = str_extract(time_slice, "[0-9]*m"),
           minute = str_remove(minute, "m"),
           minute = as.integer(minute),
           elapsed_minutes = day*24*60+hour*60+minute,
           elapsed_minutes = case_when(elapsed_minutes%%30 == 0 ~elapsed_minutes,
                                       elapsed_minutes%%30 == 3 ~elapsed_minutes-3),
           drugs =  paste(Drug1 ,Drug2, sep = "_"),
           drugs = str_remove_all(drugs, "_none|_0"),
           treatment =  paste(Drug1, Drug1Concentration,Drug2, Drug2Concentration, sep = "_"),
           treatment = str_remove_all(treatment, "_none|_0")) %>%
    filter(!elapsed_minutes %in% c(0, 30, 60, 90, 120, 150, 180, 210, 240)) %>%
    select(-matches("Ligand*|ECMp*|Endpoint*|*media")) %>%
    drop_na() %>%
    group_by(plateID, well, field, label) %>%
    mutate(migration_distance = calc_distance(.data[['Nuclei_CKn_NR_centroid-0']], .data[['Nuclei_CKn_NR_centroid-1']]),
           migration_direction = calc_direction(.data[['Nuclei_CKn_NR_centroid-0']], .data[['Nuclei_CKn_NR_centroid-1']])) %>%
    ungroup()
  
  #add in lineage tracks data
  tracks_df <- map_dfr(dir(paste0(data_path, plateID,"/Analysis/",pipeline_name,"/intermediate_files/tracking/"), pattern = "tracks.csv", recursive = TRUE, full.names = TRUE),
                       process_tracks_file)
  l1 <- l1_data %>%
    left_join(tracks_df, by = c("plateID", "well", "field", "label"))
  
  write_csv(l1,l1_data_path)
  
  set.seed(42)
  l1_subset <- l1 %>%
    group_by(well) %>%
    slice_sample(prop = .1) %>%
    ungroup()
  
  return(list(l1 = l1, l1_subset = l1_subset))
}



PCA_analysis <- function(plateID){
  l1_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_level_1.csv")
  
  df <- datasets[[plateID]][["l1_subset"]]
  pca_obj <- df %>%
    select(matches(pipeline_name), matches("neighborhood"), any_of(c("begins", "ends", "length"))) %>%
    select(-contains("centroid")) %>%
    princomp()
  scree <- screeplot(pca_obj)
  
  df_pca <- df %>%
    bind_cols(data.frame(pca_obj$scores[,1:10]))
  
  print(paste("writing level 1 subset file for ",plateID))
  write_csv(df_pca, str_replace(l1_data_path, "level_1", "level_1_subset"))
  
  df_pca_subset <- df_pca %>%
    group_by(plateID, well) %>%
    slice_sample(prop = .01) %>%
    ungroup()
  
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_PCA_plots.pdf"))
  p <- ggplot(df_pca_subset, aes(Comp.1, Comp.2, color = drugs)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(Comp.1, Comp.2, color = cell_cycle_state)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(Comp.1, Comp.2, color = length)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(Comp.1, Comp.2, color = Cell_CKn_CC_mean_intensity_ratio)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(Comp.1, Comp.2, color = lineage)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(Comp.1, Comp.2, color = migration_distance)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p_contour <- ggplot(df_pca_subset, aes(x=Comp.1, y=Comp.2) ) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
    theme_bw()
  print(p_contour)
  
  p <- ggplot(df_pca_subset, aes(x=Comp.1, y=Comp.2) ) +
    geom_density_2d_filled() +
    geom_density_2d(size = 0.25, colour = "black", bins = 25) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(level = "none",
           fill = "none") +
    theme_bw()
  print(p)
  
  res <- dev.off()
}


create_density_plots <- function(df){
  l1_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_level_1.csv")
  
  cols <- c("BEZ235_1" = "cornflowerblue", "Cabozantinib_50" = "cornflowerblue", "Trametinib_50" = "cornflowerblue",
            "BEZ235_2.5" = "darkcyan", "Cabozantinib_100" = "darkcyan","Trametinib_100" = "darkcyan",
            "BEZ235_5" = "cadetblue4",   "Cabozantinib_250" = "cadetblue4", "Trametinib_250" = "cadetblue4",
            "BEZ235_10" = "cadetblue", "Cabozantinib_500" = "cadetblue",  "Trametinib_500" = "cadetblue", 
            "BEZ235_20" = "blue",  "Cabozantinib_1000" = "blue", "Trametinib_1000" = "blue",
            "BEZ235_30" = "blue3",  "Cabozantinib_2500" = "blue3", "Trametinib_2500" = "blue3",
            "BEZ235_50" = "darkblue", "Cabozantinib_5000" = "darkblue", "Trametinib_5000" = "darkblue",
            "5FU_125" = "cornflowerblue", "5FU_250"= "darkcyan", "5FU_400" = "cadetblue4", "5FU_550" = "cadetblue",
            "5FU_750" = "blue", "5FU_1000" = "blue3", "5FU_2500" = "darkblue", 
            "AZD5438_50"= "cornflowerblue", "AZD5438_100"= "darkcyan", "AZD5438_250" = "cadetblue4",
            "AZD5438_500" = "cadetblue", "AZD5438_1000" = "blue", "AZD5438_2500" = "blue3",
            "AZD5438_5000" = "darkblue",
            "Bortezomib_1" = "cornflowerblue", "Bortezomib_2.5"= "darkcyan", "Bortezomib_5" = "cadetblue4",
            "Bortezomib_10" = "cadetblue","Bortezomib_20" = "blue", "Bortezomib_30" = "blue3","Bortezomib_50" = "darkblue",
            "Doxorubicin_0.5" = "cornflowerblue", "Doxorubicin_1"= "darkcyan", "Doxorubicin_10" = "cadetblue4",
            "Doxorubicin_25" = "cadetblue","Doxorubicin_50" = "blue", "Doxorubicin_125" = "blue3","Doxorubicin_250" = "darkblue",
            "Panobinostat.5"= "cornflowerblue", "Panobinostat_1"= "darkcyan", "Panobinostat_2.5" = "cadetblue4",
            "Panobinostat_5" = "cadetblue", "Panobinostat_6.5"  = "blue",  "Panobinostat_10" = "blue3",
            "Panobinostat_12.5" = "darkblue",
            "JQ1_50" = "cornflowerblue", "JQ1_100"= "darkcyan", "JQ1_250" = "cadetblue4", "JQ1_500" = "cadetblue",
            "JQ1_1000" = "blue", "JQ1_2500" = "blue3", "JQ1_5000" = "darkblue", 
            "MK1775_50" = "cornflowerblue", "MK1775_100"= "darkcyan", "MK1775_175" = "cadetblue4", "MK1775_275" = "cadetblue",
            "MK1775_375" = "blue", "MK1775_500" = "blue3", "MK1775_700" = "darkblue",
            "MG132_1" = "cornflowerblue", "MG132_2.5" = "darkcyan", "MG132_5" = "cadetblue4","MG132_10" = "cadetblue",
            "MG132_20" = "blue", "MG132_30" = "blue3", "MG132_50" = "darkblue",
            "Everolimus_50" = "cornflowerblue", "Everolimus_100"= "darkcyan", "Everolimus_250" = "cadetblue4",
            "Everolimus_500" = "cadetblue", "Everolimus_1000" = "blue", "Everolimus_2500" = "blue3",
            "Everolimus_5000" = "darkblue",
            "Gemcitabine_0.25" = "cornflowerblue", "Gemcitabine_0.5"= "darkcyan", "Gemcitabine_1" = "cadetblue4",
            "Gemcitabine_1.5" = "cadetblue","Gemcitabine_2"= "blue", "Gemcitabine_2.5" = "blue3","Gemcitabine_3" = "darkblue")

  #look at density plots of the cell cycle ratios per image
  #start with just the initial time point
  earliest_time_point <- unique(df$elapsed_minutes) %>%
    min() 
  time_interval <- 1440
  unique_timepoints <- unique(df$elapsed_minutes)
  selected_timepoints <- time_interval*as.integer(unique_timepoints/time_interval) %>%
    unique()
  
  df_select <- df %>%
    filter(elapsed_minutes %in% selected_timepoints) %>%
    mutate(elapsed_hours = as.integer(elapsed_minutes)/60)
  
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_density_plots.pdf"))
  
  p_densities <- ggplot(df_select, aes(Cell_CKn_CC_mean_intensity_ratio, fill = treatment, color = treatment)) +
    geom_density(alpha = .2)+
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    facet_grid(elapsed_hours~drugs) +
    guides(fill = "none", color = "none") +
    theme_bw()
  
  if("cell_cycle_state_threshold" %in% colnames(df_select)){
    p_densities <- p_densities +
      geom_vline(xintercept = cell_cycle_state_threshold)
  }
  print(p_densities)
res <- dev.off()
  }

###########

data_path <-  "/home/exacloud/gscratch/HeiserLab/images/"
pipeline_name <- "CKn"
plateIDs <- c("AU02001" = "AU02001")

datasets <- map(plateIDs, read_plate_l1)
res <- map(plateIDs, PCA_analysis)
res <- map(plateIDs, create_density_plots)

# 
# #Calculate umap object if it doesn't already exist
# #Likely want to sample the data to get a reasonable run time
# umap_obj <- umap(pca_obj$scores[,1:10])
# 
# l1_dr <- l1_subset %>%
#   bind_cols(data.frame(pca_obj$scores[,1:10]),  data.frame(umap_obj$layout))
#   
# write_csv(l1_dr, file = paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_umap.csv"))
# 
# p <- ggplot(l1_dr, aes(X1, X2)) +
#   geom_point()
# p
# 
# p <- ggplot(l1_dr, aes(X1, X2, color = slice)) +
#   geom_point()
# p
# 
# p <- ggplot(l1_dr, aes(X1, X2, color = Nuclei_CKn_NR_area )) +
#   geom_point(alpha = .5, size = .8)
# p
# 
# p <- ggplot(l1_dr, aes(X1, X2, color = begins )) +
#   geom_point(alpha = .5, size = .8)
# p
# 
# p <- ggplot(l1_dr, aes(X1, X2, color = length )) +
#   geom_point(alpha = .5, size = .8)
# p
# 
# p <- ggplot(l1_dr, aes(X1, X2, color = neighborhood_70 )) +
#   geom_point(alpha = .5, size = .8)
# p
# 
# p <- ggplot(l1_dr, aes(X1, X2, color = drugs)) +
#   geom_point(alpha = .5, size = .8)
# p

