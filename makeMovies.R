#Label and tiem stamp movies from a set of directories
library(tidyverse)
library(EBImage)
library(optparse)

#use command line arguments to identify the plate ID

#look in an Analysis directory for a csv file that holds the labels

#Get a list of the well_locations
#merge in the labels and time interval with the well_locations
#loop through all well_locations
#read in the images
#add the labels
#add the timestamps
#write the movies to the Analysis directory
