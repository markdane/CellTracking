library(tidyverse)
library(ComplexHeatmap)
library(circlize)

dataset_name = "AU565_DS2"

AU565_library <- read_csv(paste0("../metadata/",dataset_name,"_library.csv")) %>%
  group_by(plate_name, row, column, field, channel_name) %>%
  summarise(n = n(),.groups = "drop") %>%
  pivot_wider(names_from = c("column", "channel_name"), values_from = n)

AU565_library_dm <- AU565_library %>%
  select(where(is.integer)) %>%
  as.matrix()
rownames(AU565_library_dm) <- paste(AU565_library$plate_name, AU565_library$row, AU565_library$field, sep="_")

hm <- Heatmap(AU565_library_dm,
              col = colorRamp2(c(0, 190,193), c("blue", "white", "red")),
              #col = structure(2:4, names = c("193","191",   "49")), # red, green, blue
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 4),
              name = "file count")
hm
pdf(paste0("../plots/",dataset_name,"_input_image_file_count_heatmap.pdf"), height = 18)
print(hm)
result <- dev.off()

AU565_library_size <- read_csv(paste0("../metadata/",dataset_name,"_library.csv")) %>%
  group_by(plate_name, row, column, field, channel_name) %>%
  summarise(size = mean(size),.groups = "drop") %>%
  pivot_wider(names_from = c("column", "channel_name"), values_from = size)

AU565_library_size_dm <- AU565_library_size %>%
  select(matches("[12345]_*")) %>%
  as.matrix()
rownames(AU565_library_size_dm) <- paste(AU565_library_size$plate_name, AU565_library_size$row, AU565_library_size$field, sep="_")

hm_size <- Heatmap(AU565_library_size_dm,
              #col = colorRamp2(c(100, 190,193), c("blue", "white", "red")),
              #col = structure(2:4, names = c("193","192",   "49")), # red, green, blue
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 7),
              name = "average file size")
hm_size

pdf(paste0("../plots/",dataset_name,"_image_file_size_heatmaps.pdf"), height = 18)
print(hm_size)
result <- dev.off()

#Visualize combined files counts
AU565_combined <- read_csv(paste0("../metadata/",dataset_name,"_intermediate_combined_images.csv")) %>%
  group_by(plate_name, row, column) %>%
  summarise(n = n(),.groups = "drop") %>%
  pivot_wider(names_from = c("column"), values_from = n)

AU565_combined_dm <- AU565_combined %>%
  select(where(is.integer)) %>%
  as.matrix()
rownames(AU565_combined_dm) <- paste(AU565_combined$plate_name, AU565_combined$row,  sep="_")

hm_combos <- Heatmap(AU565_combined_dm,
              col = colorRamp2(c(0, 4), c("blue", "red")),
              #col = structure(2:4, names = c("193","192",   "49")), # red, green, blue
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 11),
              name = "combined file count")
hm_combos
pdf(paste0("../plots/",dataset_name,"_combined_file_count_heatmap.pdf"), height = 18)
print(hm_combos)
result <- dev.off()

#Visualize mask file counts
AU565_masks <- read_csv(paste0("../metadata/",dataset_name,"_intermediate_mask_images.csv")) %>%
  group_by(plate_name, row, column) %>%
  summarise(n = n(),.groups = "drop") %>%
  pivot_wider(names_from = c("column"), values_from = n)

AU565_masks_dm <- AU565_masks %>%
  select(where(is.integer)) %>%
  as.matrix()
rownames(AU565_masks_dm) <- paste(AU565_masks$plate_name, AU565_masks$row,  sep="_")

hm_masks <- Heatmap(AU565_masks_dm,
              col = colorRamp2(c(0, 4), c("blue", "red")),
              #col = colorRamp2(c(100, 190,193), c("blue", "white", "red")),
              #col = structure(2:4, names = c("193","192",   "49")), # red, green, blue
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 11),
              name = "mask file count")
hm_masks
pdf(paste0("../plots/",dataset_name,"_masks_file_count_heatmap.pdf"), height = 18)
print(hm_masks)
result <- dev.off()
