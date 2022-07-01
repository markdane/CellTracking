library(tidyverse)
library(rlang)
library(viridis)

calc_distance <- function(x, y){
  if(length(x) == 1 | length(y) == 1) return(d=0)
  d <- numeric(length(x)) #intialize a vector with 0s 
  for(i in 2:length(x)) { #start distances in 2nd element
    d[i] <- sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2)
  }
  return(round(d, digits = 2))
}

calc_direction <- function(x, y){
  if( length(x) == 1 | length(y) == 1) return(d=0)
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
                       plateID = col_character(),
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

create_l2_dataset <- function(l1, l1_data_path){
  l2_data_path <- str_replace(l1_data_path, "level_1", "level_2")

  l2 <- l1 %>%
    group_by(plateID, well, field, elapsed_minutes, hours, treatment, drugs, Drug1,
             Drug1Concentration, Drug2, Drug2Concentration) %>%
    summarise(n = n(),
              Nuclei_CC_mean_intensity= mean(Nuclei_CC_mean_intensity, na.rm = TRUE),
              Cyto_CC_mean_intensity= mean(Cyto_CC_mean_intensity, na.rm = TRUE),
              Cell_CC_mean_intensity_ratio= mean(Cell_CC_mean_intensity_ratio, na.rm = TRUE),
              Cell_CC_max_intensity_ratio= mean(Cell_CC_max_intensity_ratio, na.rm = TRUE),
              Cell_CC_min_intensity_ratio= mean(Cell_CC_min_intensity_ratio, na.rm = TRUE),
              G1_proportion = case_when(!n == 0  ~ sum(cell_cycle_state=="G1")/n,
                                        TRUE ~ 0),
              .groups = "drop")

  #normalize to earliest time point
  
  earliest_time_point <- unique(l2$elapsed_minutes) %>%
    min()
  
  l2 <- l2 %>%
    filter(elapsed_minutes == earliest_time_point) %>%
    rename(n_0 = n) %>%
    select("plateID","well", "field", "n_0") %>%
    right_join(l2,by = c("plateID", "well", "field")) %>%
    mutate(cell_count_norm_T0 = n/n_0)
  
  write_csv(l2,l2_data_path)
  return(l2)
}

read_plate_l1 <- function(plateID){
  l1_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_",pipeline_name,"_level_1.csv")
  if(!file.exists(l1_data_path)){
    print(paste("creating new level 1 file for ",plateID))
    l1_data <- map_dfr(dir(paste0(data_path,plateID,"/Analysis/",pipeline_name),pattern = paste0(plateID,"_.*_level_1.csv"), full.names = TRUE), read_csv, show_col_types = FALSE, lazy = FALSE,
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
                                       incubation_media = col_character(),
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
             treatment = str_remove_all(treatment, "_none|_0$|_0_none"),
             hours = elapsed_minutes/60) %>%
      select(-matches("Ligand*|ECMp*|Endpoint*|*media")) %>%
      drop_na() %>%
      group_by(plateID, well, field, label) %>%
      mutate(migration_distance = calc_distance(.data[['Cyto_CC_centroid-0']], .data[['Cyto_CC_centroid-1']]),
             migration_direction = calc_direction(.data[['Cyto_CC_centroid-0']], .data[['Cyto_CC_centroid-1']])) %>%
      ungroup()
    
    #add cell cycle state calls to each cell
    # browser()
    # untreated_ratios <- l1_data %>%
    #   filter(treatment == "Untreated") %>% 
    #   pull(Cell_CC_mean_intensity_ratio)
    # 
    # #use kmeans clustering on the Untreated cells to define a threshold for gating all cells
    # cls <- kmeans(untreated_ratios, centers = c(quantile(untreated_ratios,.25), quantile(untreated_ratios,.75)))$cluster #force lower ratio to be in cluster 1
    # cell_cycle_state_threshold <- max(untreated_ratios[cls == 1]) #threshold is the boundry between clusters
    ######Notice threshold set to 1
    cell_cycle_state_threshold <- 1
    l1_data$cell_cycle_state_threshold <- cell_cycle_state_threshold
    l1_data$cell_cycle_state <- "G1"
    l1_data$cell_cycle_state[l1_data$Cell_CC_mean_intensity_ratio>cell_cycle_state_threshold] <- "S/G2"
    
    #Add lineage values to each cell
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
    write_csv(l1_subset, str_replace(l1_data_path, "_level_1", "_level_1_subset"))
  }  else {
    l1 <- read_csv(l1_data_path, show_col_types = FALSE, lazy = FALSE,
             col_types =cols(.default = "d",
                             well = col_character(),
                             time_slice = col_character(),
                             is_parent= col_logical(),
                             plateID = col_character(),
                             Well = col_character(),
                             Cellline = col_character(),
                             Drug1 = col_character(),
                             Drug1Concentration = col_character(),
                             Drug1ConcentrationUnits = col_character(),
                             Drug2 = col_character(),
                             Drug2Concentration = col_character(),
                             Drug2ConcentrationUnits = col_character(),
                             row = col_character(),
                             treatment = col_character(),
                             drugs = col_character(),
                             cell_cycle_state = col_character()))
    l1_subset <- read_csv(str_replace(l1_data_path, "_level_1", "_level_1_subset"),
                    show_col_types = FALSE, lazy = FALSE,
                    col_types =cols(.default = "d",
                                    well = col_character(),
                                    time_slice = col_character(),
                                    is_parent= col_logical(),
                                    plateID = col_character(),
                                    Well = col_character(),
                                    Cellline = col_character(),
                                    Drug1 = col_character(),
                                    Drug1Concentration = col_character(),
                                    Drug1ConcentrationUnits = col_character(),
                                    Drug2 = col_character(),
                                    Drug2Concentration = col_character(),
                                    Drug2ConcentrationUnits = col_character(),
                                    row = col_character(),
                                    treatment = col_character(),
                                    drugs = col_character(),
                                    cell_cycle_state = col_character()))
  }
  
  l2 <- create_l2_dataset(l1, l1_data_path = l1_data_path)
  return(list(l1 = l1, l1_subset = l1_subset, l2 = l2))
}


PCA_analysis <- function(plateID){
  l1_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_level_1.csv")
  
  df <- datasets[[plateID]][["l1_subset"]]
  
  pca_obj <- df %>%
    select(matches("_CC_"), matches("neighborhood"), any_of(c("begins", "ends", "length"))) %>%
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
  
  p <- ggplot(df_pca_subset, aes(PC1, PC2, color = Cell_CC_mean_intensity_ratio)) +
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
  l1_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_level_1.csv")
  
  df <- datasets[[plateID]][["l1_subset"]]
  
  #start with just the initial time point
  earliest_time_point <- unique(df$elapsed_minutes) %>%
    min() 
  time_interval <- 1440 #one day in miinutes
  unique_timepoints <- unique(df$elapsed_minutes)
  selected_timepoints <- time_interval*as.integer(unique_timepoints/time_interval) %>%
    unique()
  
  df_select <- df %>%
    filter(elapsed_minutes %in% selected_timepoints) %>%
    mutate(hours = as.integer(elapsed_minutes)/60)
  
  df_untreated <- df_select %>%
    filter(treatment %in% c("vehicle", "control", "Untreated"))
  
  p_untreated <- ggplot(df_untreated, aes(Cell_CC_mean_intensity_ratio)) +
    geom_density(alpha = .2, color = "royalblue", fill = "royalblue")+
    theme_bw()
  
  p_untreated_by_time <- ggplot(df_untreated, aes(Cell_CC_mean_intensity_ratio, fill = treatment, color = treatment)) +
    geom_density(alpha = .2)+
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    guides(fill = "none", color = "none") +
    theme_bw()+
    facet_wrap(~hours)
  
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_",pipeline_name,"_control_density_plots.pdf"))
  print(p_untreated)
  print(p_untreated_by_time)
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
    scale_fill_manual(values=rep(viridis(8),times=400)) +
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
    select(plateID, well, field, slice, time_slice, elapsed_minutes, day,hour, minute, treatment, label, begins, ends, parent, length, lineage, Drug1, Drug1Concentration, Cell_CC_mean_intensity_ratio,cell_cycle_state) %>%
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
    scale_color_manual(values=rep(viridis(8)),times=10) +
    guides(color = "none") +
    theme_classic() +
    theme(axis.text.y = element_blank())
  p_lineages_paths

}

show_cell_cycle_plots <- function(plateID){
  df <- datasets[[plateID]][["l1"]]
  
  #dynamically assign consistent colors to the dosage concentrations
  treatments <-df$treatment%>%
    unique()  %>%
    str_sort(numeric = TRUE)
  idx <- which(treatments %in% c("control","vehicle","Untreated")) # Positions of Untreated in df$treatement
  treatments <- treatments[-idx]
  treatments_colors <-  rep(viridis(7, direction = -1), length.out = length(treatments))
  names(treatments_colors) <- treatments
  cols <- c("control"="royalblue","vehicle"="royalblue","Untreated"="royalblue", treatments_colors)
  
  df$treatment <- factor(df$treatment, levels = c(treatments, "control","vehicle","Untreated"))
  
  set.seed(42)
  df_selected<- df %>%
    filter(well %in% c("A1", "D6")) %>%
    filter(migration_distance  < 5,
           length>10) %>%
    mutate(field_label = paste0(field, "_", label)) %>%
    group_by(plateID, well) %>%
    filter(field_label %in% sample(unique(field_label), size = 5, replace = FALSE)) %>%
    group_by(plateID, well, field, label) %>%
    filter(elapsed_minutes == min(elapsed_minutes)) %>%
    select(label,plateID, well, field, 
           elapsed_minutes,
           Cyto_CC_area,
           cell_cycle_state) %>%
    mutate(t0_Cyto_CC_area = Cyto_CC_area) %>%
    select(plateID, well, field, label, t0_Cyto_CC_area) %>%
    left_join(df, by = c("plateID", "well", "field", "label"))  %>%
    mutate(Cyto_CC_area_t0norm = Cyto_CC_area/t0_Cyto_CC_area,
           label = as.character(label)) %>%
    arrange(well, field, label)
  
  p_cell_cycle_states <- ggplot(df_selected, aes(elapsed_minutes, Cell_CC_mean_intensity_ratio, color = cell_cycle_state)) +
    geom_point(aes(y = Cyto_CC_area_t0norm, group = label), color = "blueviolet") +
    geom_point(aes(y = migration_distance, group = label), color = "brown", alpha = .5) +
    geom_point(aes(group = label), size = .7) +
    labs(title = "Cell cycle measurements",
         subtitle="purple: area\nbrown: migration") +
    xlab("time") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90)) +
    coord_cartesian(ylim = c(0,10)) +
    theme_bw() +
    theme(strip.background = NULL,
          strip.text = element_text(size = 6)) +
    facet_wrap(vars(well, field, label), ncol = 5, labeller = label_wrap_gen(multi_line=FALSE))
  p_cell_cycle_states
  
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_cell_cycle_plot.pdf"), height = 30,useDingbats = FALSE)
  print(p_cell_cycle_states)
  res <- dev.off()
  
  df_lineages <- df %>%
    group_by(plateID, well, field, treatment,  label) %>%
    summarise(lifetime = mean(length-1)*.5, 
              .groups="drop")
  
  p_lineage_length <- ggplot(df_lineages, aes(x = factor(field), y = lifetime)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    xlab("field") +
    ylab("cell lifetime (hours)") +
    theme_bw()+
    facet_wrap(~treatment, ncol = 7)
  
  p_lineage_length
  
  pdf(paste0(data_path, plateID, "/Analysis/",pipeline_name,"/plots/",plateID,"_lineage_length_plot.pdf"),useDingbats = FALSE)
  print(p_lineage_length)
  res <- dev.off()
}

generate_l2_plots <- function(plateID){
  l2_data_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/",plateID,"_",pipeline_name,"_level_2.csv")
  
  df <- read_csv(l2_data_path) %>%
    mutate(field = as.factor(field)) %>% 
    arrange(elapsed_minutes)
  
  #dynamically assign consistent colors to the dosage concentrations
  treatments <-df$treatment%>%
    unique()  %>%
    str_sort(numeric = TRUE)
  idx <- which(treatments %in% c("control","vehicle","Untreated")) # Positions of Untreated in df$treatment
  treatments <- treatments[-idx]
  treatments_colors <-  rep(viridis(7, direction = -1), length.out = length(treatments))
  names(treatments_colors) <- treatments
  cols <- c("control"="royalblue","vehicle"="royalblue","Untreated"="royalblue", treatments_colors)
  
  p_n <- ggplot(df, aes(x = hours, y = n, color = field, group = field)) +
    geom_path() +
    stat_summary(fun=mean, aes(group=1), geom="line", colour="blue") +
    scale_y_continuous(trans = 'log2')+
    scale_x_continuous(breaks = c(0,24,48,72,96)) +
    labs(title = "Cell counts in each field",
         subtitle = "organized by position in the plate") +
    xlab("hours") +
    ylab("cell count (log2 spacing)") +
    facet_wrap(~well, ncol = 6)
  p_n
  
  p_n_T0 <- ggplot(df, aes(x = hours, y = cell_count_norm_T0, color = factor(field))) +
    geom_path() +
    stat_summary(fun=mean, aes(group=1), geom="line", colour="blue") +
    scale_x_continuous(breaks = c(0,24,48,72,96)) +
    scale_y_continuous(trans = 'log2') +
    labs(title = "Cell counts in each field (T0 normed)",
         subtitle = "organized by position in the plate") +
    xlab("hours") +
    ylab("cell count (log2 spacing)") +
    facet_wrap(~well, ncol = 6)
  p_n_T0
  
  p_cyto_CC_ratio <- ggplot(df, aes_string(x = "hours", y = paste0("Cell_CC_mean_intensity_ratio"), color = "factor(field)")) +
    geom_path() +
    stat_summary(fun=mean, aes(group=1), geom="line", colour="blue") +
    scale_x_continuous(breaks = c(0,24,48,72,96)) +
    labs(title = "Cell cycle intensity in the cytoplasm/nuclei, mean") +
    xlab("hours") +
    ylab("intensity") +
    facet_wrap(~well, ncol = 2)
  p_cyto_CC_ratio
  
  #summarize each well and group plots by drug
  df_control_mean <- df %>%
    filter(treatment %in% c("vehicle", "control", "Untreated")) %>%
    group_by(elapsed_minutes) %>%
    summarise(cell_count_norm_T0_control = mean(cell_count_norm_T0),
              n_control = mean(n),
              G1_proportion_control = mean(G1_proportion),
              .groups = "drop")
  
  df_mean <- df %>%
    filter(!treatment %in% c("vehicle", "control", "Untreated")) %>%
    mutate(treatment = fct_drop(treatment)) %>%
    group_by(elapsed_minutes, treatment, drugs, Drug1Concentration) %>%
    summarise(cell_count_norm_T0 = mean(cell_count_norm_T0),
              n = mean(n),
              G1_proportion = mean(G1_proportion), .groups = "drop") %>%
    left_join(df_control_mean, by = "elapsed_minutes")%>%
    mutate(hours = as.integer(elapsed_minutes)/60,
           hours = as.character(hours),
           hours = fct_inseq(hours, ordered = TRUE),
           cell_count_norm_control = n/n_control)
  
  df_mean_select <- df_mean %>%
    filter(hours %in% c(0, 24, 48, 72, 96)) %>%
    arrange(drugs, Drug1Concentration)
  
  #drug1_concentration_ranks <- unique(df_mean_select$Drug1ConcentrationRank)
  # rank_colors <-  rep(viridis(7, direction = -1), length.out = length(drug1_concentration_ranks))
  
  
  p_n_t_T0_mean <- ggplot(df_mean, aes(x = hours, y = cell_count_norm_T0, color = treatment, group = Drug1Concentration)) +
    geom_path() +
    geom_path(aes(y = cell_count_norm_T0_control), color = "firebrick4")+
    scale_x_discrete(breaks = c(0,24,48,72,96)) +
    scale_y_continuous(trans = 'log2') +
    scale_color_manual(values = cols) +
    labs(title = "Cell counts for each treatment (T0 normed)",
         subtitle = "organized by treatment and dosage") +
    xlab("hours") +
    ylab("cell count (log2 spacing)") +
    guides(color = "none") +
    facet_wrap(~drugs, ncol = 6)+
    theme(strip.text = element_text(size = 5)) +
    theme_bw()
  p_n_t_T0_mean
  
  p_G1_t_mean <- ggplot(df_mean, aes(x = hours, y = G1_proportion, color = treatment, group = Drug1Concentration)) +
    geom_path() +
    geom_path(aes(y = G1_proportion_control), color = "firebrick4")+
    scale_x_discrete(breaks = c(0,24,48,72,96)) +
    scale_color_manual(values = cols) +
    labs(title = "G1 proportion for each treatment",
         subtitle = "organized by treatment and dosage") +
    xlab("hours") +
    ylab("G1 proportion") +
    guides(color = "none") +
    facet_wrap(~drugs, ncol = 6)+
    theme(strip.text = element_text(size = 5)) +
    theme_bw()
  p_G1_t_mean
  
  concentration_colors <-  rep(viridis(7, direction = -1), length.out = 7)
  
  p_cell_count_dose_response_conc <- ggplot(df_mean_select, aes(x = Drug1Concentration, y = cell_count_norm_control, color = hours, group = hours)) +
    geom_path() +
    scale_color_manual(values = concentration_colors) +
    labs(title = "Dose response curves for cell count") +
    xlab("Drug concentration") +
    ylab("Cell count (normalized to control)") +   
    guides(color = guide_legend("Hours")) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0))+
    facet_wrap(~factor(drugs), ncol = 4, scales = "free_x")
  p_cell_count_dose_response_conc
  
  plot_path <- paste0(data_path,plateID,"/Analysis/",pipeline_name,"/plots/")
  if(!dir.exists(plot_path)) dir.create(plot_path)
  
  pdf(paste0(plot_path,plateID,"_",pipeline_name,"_QA.pdf"), width = 10, height = 14)
  print(p_n)
  print(p_n_T0)
  print(p_cyto_CC_ratio)
  res <- dev.off()
  
  pdf(paste0(plot_path, plateID,"_",pipeline_name,"_summary.pdf"), width = 10, height = 7)
  print(p_n_t_T0_mean)
  print(p_G1_t_mean)
  #print(p_cell_count_dose_response)
  print(p_cell_count_dose_response_conc)
  res <- dev.off()
}

###########

data_path <-  "/home/exacloud/gscratch/HeiserLab/images/"
pipeline_name <- "CtcK"
plateIDs <- c("HC01301" = "HC01301")

datasets <- map(plateIDs, read_plate_l1)
res <- map(plateIDs, generate_l2_plots)

#####These all need to be debugged
#res <- map(plateIDs, create_control_density_plots)

# res <- map(plateIDs, create_lineage_pdf, wll = "A1", fld = 2)
# res <- map(plateIDs, create_lineage_pdf, wll = "B1", fld = 2)
# res <- map(plateIDs, create_lineage_pdf, wll = "C1", fld = 2)
# res <- map(plateIDs, create_lineage_pdf, wll = "D2", fld = 2)
# res <- map(plateIDs, create_lineage_pdf, wll = "B3", fld = 2)
# res <- map(plateIDs, create_lineage_pdf, wll = "C3", fld = 2)
# res <- map(plateIDs, create_lineage_pdf, wll = "D4", fld = 2)
# res <- map(plateIDs, create_lineage_pdf, wll = "B5", fld = 2)
# res <- map(plateIDs, create_lineage_pdf, wll = "C5", fld = 2)
# res <- map(plateIDs, create_lineage_pdf, wll = "D6", fld = 2)

#res <- map(plateIDs, show_cell_cycle_plots)

