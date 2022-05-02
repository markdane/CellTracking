library(tidyverse)
library(viridis)

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
  #latch cell cycle state after 2 consecutive G2/S calls
  # determine_ccs_latch <-function(em, ratio, ccs){
  #   browser()
  #   return(ccs)
  # }
  
  # foo <- l1_data %>%
  #   group_by(plateID, well, field, label) %>%
  #   arrange(elapsed_minutes) %>%
  #   mutate(cell_cycle_state_latch = determine_ccs_latch(elapsed_minutes, Cell_CKn_CC_mean_intensity_ratio,cell_cycle_state))
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
  
  df_pca_scores <- as_tibble(pca_obj$scores) %>%
    rename_with(~gsub("Comp.", "PC", .x, fixed = TRUE)) %>%
    select(num_range("PC", 1:10))
  
  df_pca <- bind_cols(df, df_pca_scores)
  
  scree <- screeplot(pca_obj)
  
  df_pca_subset <- df_pca %>%
    group_by(plateID, well) %>%
    slice_sample(prop = .2) %>%
    ungroup()
  
  plots_path <- paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/")
  if(!dir.exists(plots_path)) dir.create(plots_path)
  
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_PCA_plots.pdf"),useDingbats = FALSE)
  
  print(screeplot(pca_obj))
  
  p <- ggplot(df_pca_subset, aes(PC1, PC2, color = drugs)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(PC1, PC2, color = cell_cycle_state)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(PC1, PC2, color = length)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(PC1, PC2, color = Cell_CKn_CC_mean_intensity_ratio)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(PC1, PC2, color = lineage)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p <- ggplot(df_pca_subset, aes(PC1, PC2, color = migration_distance)) +
    geom_point(size = .5, alpha = .3) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw()
  print(p)
  
  p_contour <- ggplot(df_pca_subset, aes(x=PC1, y=PC2) ) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
    theme_bw()
  print(p_contour)
  
  p <- ggplot(df_pca_subset, aes(x=PC1, y=PC2) ) +
    geom_density_2d_filled() +
    geom_density_2d(size = 0.25, colour = "black", bins = 25) +
    #coord_cartesian(xlim = c(-700,500), ylim = c(-300,300)) +
    guides(level = "none",
           fill = "none") +
    theme_bw()
  print(p)
  
  res <- dev.off()
}

create_density_plots <- function(plateID){
  #browser()
  l1_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_level_1.csv")
  
  df <- datasets[[plateID]][["l1_subset"]]
  

  #dynamically assign consistent colors to the dosage concentrations
  treatments <-df$treatment%>%
    unique()  %>%
    str_sort(numeric = TRUE)
  idx <- which(treatments == "Untreated") # Positions of Untreated in df$treatement
  treatments <- treatments[-idx]
  treatments_colors <-  rep(viridis(7, direction = -1), length.out = length(treatments))
  names(treatments_colors) <- treatments
  cols <- c("Untreated"="royalblue", treatments_colors)
  
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
  
  plateID <- unique(df_select$plateID)
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_density_plots.pdf"))
  
  p_densities <- ggplot(df_select, aes(Cell_CKn_CC_mean_intensity_ratio, fill = treatment, color = treatment)) +
    geom_density(alpha = .2)+
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    facet_grid(elapsed_hours~drugs) +
    guides(fill = "none", color = "none") +
    theme_bw()
  
  if("cell_cycle_state_threshold" %in% colnames(df_select)){
    cell_cycle_state_threshold <- unique(df_select$cell_cycle_state_threshold)
    p_densities <- p_densities +
      geom_vline(xintercept = cell_cycle_state_threshold)
  }
  print(p_densities)
  res <- dev.off()
}


create_lineage_pdf <- function(plateID, wll, fld){
  df <- datasets[[plateID]][["l1"]]
  
  df_plateID_well_field <- df %>%
    filter(plateID == plateID,
           well ==  wll,
           field == fld)
  
  t0 <- unique(df_plateID_well_field$elapsed_minutes) %>% min()
  
  t0_lineages <- df_plateID_well_field %>%
    filter(elapsed_minutes == t0) %>%
    select(lineage) %>%
    unique() %>%
    drop_na() %>%
    #slice_sample(prop = .5) %>%
    pull()
  
  #make a complete lineage tree diagram
  df_lineage_labels <- df_plateID_well_field %>%
    filter(lineage %in%t0_lineages) %>%
    arrange(lineage) %>%
    mutate(lineage = factor(lineage),
           lineage = fct_inorder(lineage),
           lineage_label = paste0(lineage,"_",label))
  
  p_lineages_state <- ggplot(df_lineage_labels, aes(x = elapsed_minutes, y = lineage_label, group = label, fill = lineage, color = cell_cycle_state)) +
    geom_point(size = 1, shape = 21, stroke = .1)+
    scale_fill_manual(values=rep(brewer.pal(8,"Dark2"),times=400)) +
    scale_color_manual(values = c("G1" = "transparent", "S/G2" = "black"))+
    guides(color = "none", fill = "none") +
    labs(title = paste("T0 lineages -",unique(df_lineage_labels$Drug1),"field",unique(df_lineage_labels$field)),
         subtitle = "S/G2 shown with black outline") +
    ylab("Cells colored by lineage") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank()) 
  p_lineages_state
  
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_",wll,"_", fld,"_lineage_plot.pdf"), height = 30,useDingbats = FALSE)
  print(p_lineages_state)
  res <- dev.off()
}


generate_lineage_plots <- function(plateID){
#Use the l1 data to create lineage plots of each movie
  selected_field <- 1

  df <- datasets[[plateID]][["l1"]]

  t0 <- unique(df$elapsed_minutes) %>% min()

  t0_lineages <- df %>%
    filter(elapsed_minutes == t0) %>%
    select(lineage) %>%
    unique() %>%
    drop_na() %>%
    slice_sample(prop = .1) %>%
    pull()

  df_t0_cells <- df %>%
    select(plateID, well, field, slice, time_slice, elapsed_minutes, day,hour, minute, treatment, label, begins, ends, parent, length, lineage, Drug1, Drug1Concentration, Cell_CKn_CC_mean_intensity_ratio,cell_cycle_state) %>%
    filter(field == selected_field,
           label %in% t0_lineages) %>%
    arrange(length) %>%
    mutate(label = factor(label),
           label = fct_inorder(label))

  #View the generation 0 cells colored by cell state
  p_gen_0_cells <- ggplot(df_t0_cells, aes(x = elapsed_minutes, y = label, color = cell_cycle_state)) +
    geom_path(size = .5) +
    scale_color_manual(values = c("G1" = "black", "S/G2" = "red")) +
    labs(title = "Generation 0 cells",
         subtitle = "colored by cell cycle state") +
    xlab("minutes") +
    ylab("") +
    theme_bw() +
    theme(axis.ticks.y=element_blank(),
          axis.text.y=element_blank())
  p_gen_0_cells

  #make a complete lineage tree diagram
  df_lineage_labels <- df %>%
    filter(lineage %in%t0_lineages) %>%
    arrange(lineage) %>%
    mutate(lineage = factor(lineage),
           lineage = fct_inorder(lineage),
           lineage_label = paste0(lineage,"_",label))

  track_lengths<- df %>%
    group_by(well, field, label) %>%
    summarise(track_length = unique(length))

  mean_track_length <- track_lengths %>%
    group_by(well, field) %>%
    summarise(mean_track_length = mean(track_length))
  # %>%
  #   filter(!(begins <3 & length < 6))

  p_lineages_paths <- ggplot(df_lineage_labels, aes(x = elapsed_minutes, y = lineage_label, group = factor(label), color = factor(lineage))) +
    geom_path() +
    scale_color_manual(values=rep(brewer.pal(8,"Dark2"),times=10)) +
    guides(color = "none") +
    theme_classic() +
    theme(axis.text.y = element_blank())
  p_lineages_paths

}

###########

data_path <-  "/home/exacloud/gscratch/HeiserLab/images/"
pipeline_name <- "CKn"
plateIDs <- c("AU03501" = "AU03501")

datasets <- map(plateIDs, read_plate_l1)
res <- map(plateIDs, PCA_analysis)
res <- map(plateIDs, create_density_plots)

res <- map(plateIDs, create_lineage_pdf, wll = "A1", fld = 2)
res <- map(plateIDs, create_lineage_pdf, wll = "B1", fld = 2)
res <- map(plateIDs, create_lineage_pdf, wll = "C1", fld = 2)
res <- map(plateIDs, create_lineage_pdf, wll = "D2", fld = 2)
res <- map(plateIDs, create_lineage_pdf, wll = "B3", fld = 2)
res <- map(plateIDs, create_lineage_pdf, wll = "C3", fld = 2)
res <- map(plateIDs, create_lineage_pdf, wll = "D4", fld = 2)
res <- map(plateIDs, create_lineage_pdf, wll = "B5", fld = 2)
res <- map(plateIDs, create_lineage_pdf, wll = "C5", fld = 2)
res <- map(plateIDs, create_lineage_pdf, wll = "D6", fld = 2)

show_cell_cycle_plots <- function(plateID){
  df <- datasets[[plateID]][["l1"]]
  
  set.seed(42)
  df_selected<- df %>%
    group_by(plateID, well, field) %>%
    filter(label %in% sample(unique(label), size = 5, replace = FALSE)) %>%
    group_by(plateID, well, field, label) %>%
    filter(elapsed_minutes == min(elapsed_minutes)) %>%
    select(label,
           elapsed_minutes,
           Nuclei_CKn_NR_area,
           cell_cycle_state) %>%
    mutate(t0_Nuclei_CKn_NR_area = Nuclei_CKn_NR_area) %>%
    select(plateID, well, field, label, t0_Nuclei_CKn_NR_area) %>%
    left_join(df, by = c("plateID", "well", "field", "label"))  %>%
    filter(migration_distance  < 5,
           length>10) %>%
    mutate(Nuclei_CKn_NR_area_t0norm = Nuclei_CKn_NR_area/t0_Nuclei_CKn_NR_area)
  
  p_cell_cycle_states <- ggplot(df_selected, aes(elapsed_minutes, Cell_CKn_CC_mean_intensity_ratio, color = cell_cycle_state)) +
    geom_path(aes(y = Nuclei_CKn_NR_area_t0norm), color = "blueviolet") +
    geom_path(aes(y = migration_distance), color = "brown", alpha = .5) +
    geom_point(size = 1) +
    labs(title = "Cell cycle measurements",
         subtitle="purple: area\nbrown: migration") +
    xlab("time") +
    ylab("") +
    facet_wrap(~label)
  
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_cell_cycle_plot.pdf"), height = 30,useDingbats = FALSE)
  print(p_cell_cycle_states)
  res <- dev.off()
}
res <- map(plateIDs, show_cell_cycle_plots)

