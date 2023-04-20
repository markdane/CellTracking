# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# R-Script for animating time-series data
# Based on a script from Joachim Goedhart (@joachimgoedhart)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(gganimate)
library(magick)
library(readxl)
library(gridExtra)

#functions####

#read and prepare data
read_process_data <- function(combo_name){
  drug_names <- str_split(combo_name, pattern = "_", simplify = TRUE) %>%
    tolower()
  
  cell <- read_excel(paste0("../data/",combo_name,"_CELL.xlsx")) %>%
    janitor::clean_names() %>%
    pivot_longer(cols = !time, names_to = "treatment", values_to = "counts")
  
  cycle <- read_excel(paste0("../data/",combo_name,"_CYCLE.xlsx")) %>%
    janitor::clean_names()%>%
    pivot_longer(cols = !time, names_to = "treatment", values_to = "G1_percent")
  
  return(list(cell = cell, cycle = cycle))
}

#Create static plots to help development
plot_data <- function(cell_cycle_data, exclude_combo = TRUE){
  nframes <- length(unique(cell_cycle_data$cell$time))
  
  #Choose a dark palette, and then replace the control and combo colors
  #replace dark gray with black
  treatment_colors <- RColorBrewer::brewer.pal(name = "Dark2", n = length(unique(cell_cycle_data$cell$treatment)))
  names(treatment_colors) <- unique(cell_cycle_data$cell$treatment)
  treatment_colors[["control"]] <- "black"
  treatment_colors[["combo"]] <- "blue"
  
  if(exclude_combo){
    df <- cell_cycle_data$cell %>%
      filter(!treatment == "combo")
  } else {
    df <- cell_cycle_data$cell
  }
  
  # Code for an ordinary plot
  # p <- ggplot(df, aes(x=time, y=counts, color=treatment)) + 
  #   geom_line(size=.5) + 
  #   scale_x_continuous(breaks = c(0,24,48, 72, 96)) +
  #   scale_color_manual(values = treatment_colors) +
  #   lims(y = c(0,3.5)) +
  #   labs(x = "Time (hours)",
  #        y = "Cell counts (fold change)",
  #        color = "Treatments") +
  #   theme_light(base_size = 16)+
  #   theme(panel.background = element_blank(),
  #         #panel.border = element_blank(),
  #         #panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()
  #   )
  # p
  
  # Generate the cell count animation
  if(exclude_combo) {
    cell_plot <- ggplot(df, aes(x=time, y=counts, color=treatment)) + 
      geom_line(size=.7) +
      geom_point(size=1) +
      transition_reveal(time) +
      scale_x_continuous(breaks = c(0,24,48, 72, 96)) +
      scale_color_manual(values = treatment_colors) +
      lims(y = c(0,3.5)) +
      labs(x = "Time (hours)",
           y = "Cell counts (fold change)",
           color = "Treatments") +
      theme_light(base_size = 16)+
      theme(panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="none"
      )
  } else {
    df_no_combo <- df %>%
      filter(!treatment == "combo")
    df_combo <- df  %>%
      filter(treatment == "combo")
    p <- ggplot(df_no_combo, aes(x=time, y=counts, color=treatment)) + 
      geom_line(size=.7, alpha = .5) +
      shadow_mark(past = TRUE, future = TRUE, exclude_layer = c(2,5)) +
      scale_x_continuous(breaks = c(0,24,48, 72, 96)) +
      scale_color_manual(values = treatment_colors) +
      lims(y = c(0,3.5)) +
      labs(x = "Time (hours)",
           y = "Cell counts (fold change)",
           color = "Treatments") +
      theme_light(base_size = 16)+
      theme(panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="none"
      )
    cell_plot <- p + geom_line(data = df_combo, size=2, alpha = 1) +
      transition_reveal(time) 
    
  }
  
  #start cell cycle plots
  # Code for an ordinary plot
  if(exclude_combo){
    df <- cell_cycle_data$cycle %>%
      filter(!treatment == "combo")
  } else {
    df <- cell_cycle_data$cycle
  }
  
  # p <- ggplot(df, aes(x=time, y=G1_percent, color=treatment)) +
  #   geom_line(size=.5) +
  #   scale_x_continuous(breaks = c(0,24,48, 72, 96)) +
  #   scale_color_manual(values = treatment_colors) +
  #   lims(y = c(0,100)) +
  #   labs(x = "Time (hours)",
  #        y = "G1 percent",
  #        color = "Treatments") +
  #   theme_light(base_size = 16)+
  #   theme(panel.background = element_blank(),
  #         panel.grid.minor = element_blank()
  #   )
  # p
  
  # Generate the cell cycle phase animation
  if(exclude_combo) {
    cycle_plot <- ggplot(df, aes(x=time, y=G1_percent, color=treatment)) +
      geom_line(size=.7) +
      geom_point(size=1) +
      transition_reveal(time) +
      scale_x_continuous(breaks = c(0,24,48, 72, 96)) +
      scale_color_manual(values = treatment_colors) +
      lims(y = c(0,100)) +
      labs(x = "Time (hours)",
           y = "G1 percent",
           color = "Treatments") +
      theme_light(base_size = 16)+
      theme(panel.background = element_blank(),
            panel.grid.minor = element_blank()
      )
  } else {
    df_no_combo <- df %>%
      filter(!treatment == "combo")
    df_combo <- df  %>%
      filter(treatment == "combo")
    p <- ggplot(df_no_combo, aes(x=time, y=G1_percent, color=treatment)) +
      geom_line(size=.7, alpha = .5) +
      shadow_mark(past = TRUE, future = TRUE, exclude_layer = 2) +
      scale_x_continuous(breaks = c(0,24,48, 72, 96)) +
      scale_color_manual(values = treatment_colors) +
      lims(y = c(0,100)) +
      labs(x = "Time (hours)",
           y = "G1 percent",
           color = "Treatments") +
      theme_light(base_size = 16)+
      theme(panel.background = element_blank(),
            panel.grid.minor = element_blank()
      )
    cycle_plot <- p + geom_line(data = df_combo, size=2, alpha = 1) +
      transition_reveal(time)
    
  }
  
  return(list(cell = cell_plot, cycle = cycle_plot))
}

#animate the plots
animate_plots <- function(cell_cycle_plot, combo_name){
  nframes <- length(unique(cell_cycle_plot$cell$data$time))
  panel_a <- animate(cell_cycle_plot$cell, nframes, height = 500, width = 500, renderer = magick_renderer())
  panel_b <- animate(cell_cycle_plot$cycle, nframes, height = 500, width = 610, renderer = magick_renderer())
  
  #Create first frame then add the rest
  full_gif <- image_append(c(panel_a[1], panel_b[1]))
  for(i in 2:nframes){
    next_frame <- image_append(c(panel_a[i], panel_b[i]))
    full_gif <- c(full_gif, next_frame)
  }
  #Save the combined GIF
  if(cell_cycle_plot$cell$plot_env[["exclude_combo"]]){
    image_write_gif(full_gif, paste0("../plots/",combo_name,"_no_combo.gif"))
  } else {
    image_write_gif(full_gif, paste0("../plots/",combo_name,"_dynamic_combo.gif"))
  }
}

#end of functions###

#read the data files
combo_names <- c("LAP100_PALBO100" = "LAP100_PALBO100",
                 "GEM10_PALBO100" = "GEM10_PALBO100",
                 "GEM10_LAP100" = "GEM10_LAP100",
                 "GEM10_DOX20" = "GEM10_DOX20")
cell_cycle_data <- map(combo_names, read_process_data)

#generate the plots
cell_cycle_plots <- map(cell_cycle_data, plot_data, exclude_combo = FALSE)

#animate the plots
res <- imap(cell_cycle_plots, animate_plots)

#Do a version with the combinations
cell_cycle_plots <- map(cell_cycle_data, plot_data, exclude_combo = TRUE)

#animate the plots
res <- imap(cell_cycle_plots, animate_plots)


