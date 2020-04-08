library(tidyverse)
library(here)
library(rmarkdown)

#assume there is a subfolder named Data that contains data files
#with plateID_well alphanumeric_field_intensities.csv files

process_files <- function(path){
  out_filename <- str_replace(path,"Points","intensities")
  
  if(!dir.exists("Reports")) dir.create("Reports")
  
  report_name <- str_replace(path, "Data","Reports") %>%
    str_replace("_Points.csv",".html")
  
  render("MTrackJ/Manual_tracking_EDA.Rmd",
         output_file = report_name,
         output_format = "html_document",
         output_options = list(code_folding = "hide"))
}
 
res <- dir("MTrackJ/Data", pattern = "Points.csv", full.names = TRUE)[1] %>%
  here() %>%
  map(process_files)
