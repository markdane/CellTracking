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

selected_wells <- AU01101_02_CELL %>%
  filter(plate == "AU01101",
         well %in% c("a1", "a2", "a3", "a4", "a4", "a6"))

#Read the number of frames
nframes <- length(unique(selected_wells$time))

# Code for an ordinary plot - uncomment to have to inspect the data
p <- ggplot(selected_wells, aes(x=time, y=value, color=well)) + 
  geom_line(size=.5) + 
  #geom_point(size=3) +
  scale_x_continuous(breaks = c(0,24,48, 72, 96)) +
  facet_wrap(~plate, nrow = 3, strip.position = "right") +
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        strip.background =element_rect(fill=NA, color='grey'))
p


# Generate the animation
animated_plot <- ggplot(selected_wells, aes(x=time, y=value, color=well)) + 
  geom_line(size=.5) +
  geom_point(size=1) +
  transition_reveal(time) +
  scale_x_continuous(breaks = c(0,24,48, 72, 96)) +
  facet_wrap(~plate, nrow = 3, strip.position = "right") +
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        strip.background =element_rect(fill=NA, color='grey'))

nframes <- 97
# Uncomment the next line to save the animated plot as GIF
#image_write_gif(animate(animated_plot, nframes, renderer = magick_renderer()), 'data_animation.gif')

panel_b <- animate(animated_plot, nframes, height = 1024, renderer = magick_renderer())

panel_a <- image_read("LI802303_A02_2_treatment.gif")

last_frame <- image_append(c(panel_a[nframes], panel_b[nframes]))

#debug 
nframes = 5
#Create first frame
two_movies <- image_append(c(panel_a[1], panel_a[1])) 
combined_gif <- image_append(c(two_movies, panel_b[1]), stack = TRUE)
for(i in 2:nframes){
  two_movies <- image_append(c(panel_a[i], panel_a[i])) 
  combined_gif <- image_append(c(two_movies, panel_b[i]), stack = TRUE)
}

#Save the combined GIF
image_write_gif(combined_gif, 'Combined.gif')
