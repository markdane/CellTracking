#format level 1 data  for napari tracks layer
library(tidyverse)

#start hard coded then refactor to generalize
#read in data file
data_path <- "../../AU565/images/"
plateID <- "AU02001"
well <- "A1"
field <- 1

l1 <- read_csv(paste0(data_path,plateID,"/Analysis/PC/",plateID,"_",well,"_field_",field,"_level_1.csv"))

# data (array (N, D+1)) – Coordinates for N points in D+1 dimensions. ID,T,(Z),Y,X. 
# The first axis is the integer ID of the track. D is either 3 or 4 for planar or volumetric timeseries respectively.

tracks <- l1 %>%
  select(label, slice, `Nuclei_PC_NR_centroid-1`, `Nuclei_PC_CC_centroid-0`)

# graph (dict {int: list}) – Graph representing associations between tracks. 
# Dictionary defines the mapping between a track ID and the parents of the track. 
# This can be one (the track has one parent, and the parent has >=1 child) in the case of track splitting, 
# or more than one (the track has multiple parents, but only one child) in the case of track merging. 
# See examples/tracks_3d_with_graph.py

graph <- paste0(data_path,plateID,"/Analysis/PC/intermediate_files/tracking/",well,"/field_",field,"/results/tracks.csv") %>%
  read_csv() %>%
  select(label, parent)

#properties (dict {str: array (N,)}, DataFrame) – Properties for each point. 
#Each property should be an array of length N, where N is the number of points.

properties <- l1 %>%
  select()