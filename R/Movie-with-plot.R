# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Movie with a plot
# R-Script for plotting time-series data, retrieved from FIJI and adding it to a movie
# Created by Joachim Goedhart (@joachimgoedhart), first version 2020
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(gganimate)
library(magick)
library(readxl)

#read the data file and reformat into two 24 well plates with time values from 0 to 96
AU01101_02_CELL <- read_excel("data/AU01101-02_CELL.xlsx") %>%
  janitor::clean_names() %>%
  select(-time_26) %>%
  pivot_longer(cols = !time_1) %>%
  mutate(index = str_remove(name, ".*_"),
         index = as.integer(index),
         plate = case_when(index <=25 ~ "AU01101",
                           index > 25  ~"AU01102"),
         well = str_remove(name, "_.*")) %>%
  select(-index) %>%
  rename(time = time_1)

#Read the number of frames
nframes <- length(unique(AU01101_02_CELL$time))

# Code for an ordinary plot - uncomment to have to inspect the data
p <- ggplot(AU01101_02_CELL, aes(x=time, y=value, color=well)) + 
  geom_line(size=2) + 
  geom_point(size=3) +
  facet_wrap(~plate)
p


# Generate the animation
animated_plot <- ggplot(AU01101_02_CELL, aes(x=time, y=value, color=well)) + geom_line(size=.5) + geom_point(size=1) + transition_reveal(time) +
  facet_wrap(~plate)

#Do some formatting of the layout
animated_plot <- animated_plot + theme_light(base_size = 16)
animated_plot <- animated_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
animated_plot <- animated_plot + theme(legend.position="none")
animated_plot <- animated_plot + facet_wrap(~plate, nrow = 3, strip.position = "right") + 
  theme(strip.background =element_rect(fill=NA, color='grey')) + 
  theme(strip.text = element_text(colour = 'black'))

# Uncommen the next line to save the animated plot as GIF
image_write_gif(animate(animated_plot, nframes, renderer = magick_renderer()), 'animation.gif')

panel_b <- animate(animated_plot, nframes, height = 529, renderer = magick_renderer())

panel_a <- image_read("LI802303_A02_2_treatment.gif")

combined_gif <- image_append(c(panel_a[1], panel_b[1]))

#debug 
nframes <- 97
for(i in 2:nframes){
  combined_panel <- image_append(c(panel_a[i], panel_b[i]))
  combined_gif <- c(combined_gif, combined_panel)
}

#Save the combined GIF
image_write_gif(combined_gif, 'Combined.gif')
