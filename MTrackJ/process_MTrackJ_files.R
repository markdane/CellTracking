library(tidyverse)
library(knitr)
library(rmarkdown)

#assume there is a subfolder named Data that contains data files
#with plateID_well alphanumeric_field_intensities.csv files

process_files <- function(path){
  filename <- str_remove(path, ".*/")
  plateID <-  str_remove(filename, ".*/") %>%
    str_remove( "_.*")
  well <- str_extract(filename, "_[[:alpha:]][[:digit:]]*_") %>%
    str_remove_all( "_")
  field <- str_extract(filename, "_[[:digit:]]*_") %>%
    str_remove_all("_")
  out_filename <- str_replace(path,"Points","intensities")
  if(!dir.exists("Reports")) dir.create("Reports")
  report_name <- paste0("Reports/",plateID,"_",well, "_", field,".html")
  
  render("Manual_tracking_EDA.Rmd",
         output_file = report_name,
         output_format = "html_document")
  
}
 
files_to_process <- dir("Data", pattern = "Points.csv", full.names = TRUE)

res <- map(files_to_process, process_files)
