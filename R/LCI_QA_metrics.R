#calculate and save tracking metrics

library(tidyverse)
calc_QA_metrics <- function(fn){
  df <- read_csv(fn, show_col_types = FALSE)
  unique_cells <- unique(df$label)
  unique_cells_count <- length(unique_cells)
  #print(paste("Cell count:",unique_cells_count))
  starting_cells <- unique(df$label[df$begins == 0])
  starting_cells_count <- length(starting_cells)
  #print(paste("Cell count at experiment start (founders):",starting_cells_count))
  
  starting_cells_lengths <- df %>%
    filter(label %in% starting_cells) %>%
    select(label, length) %>%
    distinct()
  #print("Founder cells' lifetimes summary (hrs)")
  
  #summary((starting_cells_lengths$length-1)*.5)
  starting_cells_lived_full_exp <- unique(df$label[df$begins == 0 & df$length == 193])
  starting_cells_lived_full_exp_count <- length(starting_cells_lived_full_exp)
  #print(paste("Number of founder cells that lived the full experiment:",starting_cells_lived_full_exp_count))
  starting_cells_parents <- unique(df$label[df$begins == 0 & df$is_parent])
  starting_cells_parents_count <- length(starting_cells_parents)
  #print(paste("Number of founder cells that divided:",starting_cells_parents_count))
  starting_cells_no_offspring <-  unique(df$label[df$begins == 0 & !df$is_parent & df$length < 193])
  starting_cells_no_offspring_count <-  length(starting_cells_no_offspring)
  #print(paste("Number of founder cells that did not divide",starting_cells_no_offspring_count))
  
  no_parent_cells <- unique(df$label[df$begins>0 & df$label != df$parent])
  no_parent_cells_count <- length(no_parent_cells)
  #print(paste("Number of orphan cells",no_parent_cells_count))
  
  no_parent_lengths <- df %>%
    filter(label %in% no_parent_cells) %>%
    select(label, length) %>%
    distinct()
  #print("Summary of lifetimes (hrs) of orphan cells")
  #summary((no_parent_lengths$length-1)*.5)
  
  if("root" %in% colnames(df)){
    lineages <- df %>%
      select(label, root) %>%
      distinct() %>%
      group_by(root) %>%
      mutate(cells_in_lineage = n()) %>%
      select(!label) %>%
      distinct()
    #print("Summary of number of cells in each lineage")
    #summary(lineages$cells_in_lineage)
    #table(lineages$cells_in_lineage, dnn = c("Count of lineages with labeled cell number"))
  }
  results <- list(filename = fn,
                  starting_cells_count = starting_cells_count,
                  starting_cells_lived_full_exp_count = starting_cells_lived_full_exp_count,
                  starting_cells_parents_count = starting_cells_parents_count,
                  starting_cells_no_offspring_count = starting_cells_no_offspring_count,
                  no_parent_cells_count = no_parent_cells_count
                  
  )
  return(results)
}

foo <- "/home/"
data_path <- "/home/groups/heiserlab_genomics/home/dane/CellTracking/images/"
plate_id <- "HC01501"
filenames <- dir(paste0(data_path,plate_id),pattern = "level",recursive = TRUE,full.names = TRUE)
results_df <- map_df(filenames, calc_QA_metrics)
res <- write_csv(results_df, file = paste0(data_path, plate_id,"/QA/QA_metrics.csv"))
#pipeline_plate_well_file_path <- "/CntB/cell_level_data/HC01501_B2_field_1_level_0.csv"
#pipeline_plate_well_file_path <- "HC00801/CtncT/cell_level_data/HC00801_A1_field_1_level_1.csv"
#pipeline_plate_well_file_path <- "HC00801/CtcB/cell_level_data/HC00801_A1_field_1_level_1.csv"
#df <- read_csv(paste0(data_path, pipeline_plate_well_file_path), show_col_types = FALSE)
