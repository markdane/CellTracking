library(tidyverse)

pipeline_name = "PI"
data_path <- "../../../../images/"
dataset_name = "AU565_DS2"
#plateIDs <- c("AU00601","AU00602","AU00701","AU00702","AU00801","AU00802","AU00901","AU00902","AU01001","AU01002","AU01101","AU01102")
plateIDs <- c("AU01401","AU01501","AU01502","AU01601","AU01602","AU01701","AU01702","AU01801","AU01802","AU01901","AU01902","AU02001","AU02002","AU02101")
plateIDs <- c("2101201","2101301")
data_paths <- paste0(data_path, plateIDs)


#read file names into a dataframe
df <- map(data_paths, dir, pattern = "m.tif", recursive = TRUE) %>%
  unlist() %>%
  bind_cols() %>%
  rename("filename" = "...1") %>%
  mutate(filename = str_remove(filename, ".*/"),
         plate_name = str_remove_all(filename, "_.*"),
         channel_name = str_extract(filename, "_[RPG]_"),
         channel_name = str_remove_all(channel_name,"_"),
         well = str_extract(filename, "_[ABCD123456]*_"),
         well = str_remove_all(well,"_"),
         row = str_extract(well, "[[:alpha:]]"),
         column = str_extract(well, "[[:digit:]*]"),
         well = paste0(row,column),
         field = str_extract(filename, "_[1234]_"),
         field = str_remove_all(field,"_"),
         size = file.size(paste0(data_path,plate_name,"/",well,"/field_",field,"/",filename)),
         timepoint = str_extract(filename, "[[:alnum:]]*.tif"),
         timepoint = str_remove_all(timepoint,".tif")
  )

write_csv(df, paste0("../metadata/",dataset_name,"_library.csv"))

#Get the number of combined image files ready for ilastik
combined_files_df <- map(data_paths, dir, pattern = "_stack.tif", recursive = TRUE) %>%
  unlist() %>%
  bind_cols() %>%
  rename("filename" = "...1") %>%
  mutate(filename = str_remove(filename, ".*/"),
         plate_name = str_remove_all(filename, "_.*"),
         channel_name = str_extract(filename, "_[RPG]*_"),
         channel_name = str_remove_all(channel_name,"_"),
         well = str_extract(filename, "_[ABCD123456]*_"),
         well = str_remove_all(well,"_"),
         row = str_extract(well, "[[:alpha:]]"),
         column = str_extract(well, "[[:digit:]*]"),
         well = paste0(row,column),
         field = str_extract(filename, "_[1234]_"),
         field = str_remove_all(field,"_"),
         size = file.size(paste0(data_path,plate_name,"/Analysis/",pipeline_name,"/intermediate_files/",filename))
)

write_csv(combined_files_df, paste0("../metadata/",dataset_name,"_intermediate_combined_images.csv"))

#Get the number of pixel mask files processed by ilastik
mask_files_df <- map(data_paths, dir, pattern = " Segmentation.h5", recursive = TRUE) %>%
  unlist() %>%
  bind_cols() %>%
  rename("filename" = "...1") %>%
  mutate(filename = str_remove(filename, ".*/"),
         plate_name = str_remove_all(filename, "_.*"),
         channel_name = str_extract(filename, "_[RPG]*_"),
         channel_name = str_remove_all(channel_name,"_"),
         well = str_extract(filename, "_[ABCD123456]*_"),
         well = str_remove_all(well,"_"),
         row = str_extract(well, "[[:alpha:]]"),
         column = str_extract(well, "[[:digit:]*]"),
         well = paste0(row,column),
         field = str_extract(filename, "_[1234]_"),
         field = str_remove_all(field,"_"),
         size = file.size(paste0(data_path,plate_name,"/Analysis/",pipeline_name,"/intermediate_files/",filename))
  )

write_csv(mask_files_df, paste0("../metadata/",dataset_name,"_intermediate_mask_images.csv"))
