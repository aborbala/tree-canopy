# Install and load necessary packages
if (!requireNamespace("sf", quietly = TRUE)) {
  install.packages("sf")
}
if (!requireNamespace("lidR", quietly = TRUE)) {
  install.packages("lidR")
}
if (!requireNamespace("tools", quietly = TRUE)) {
  install.packages("tools")
}

library(sf)
library(lidR)
library(tools)

# Function to remove buildings from Lidar data
remove_buildings_from_las <- function(buildings_path, las_path) {
  buildings <- st_read(buildings_path)
  crs <- st_crs(buildings)
  las <- readLAS(las_path, filter = "-drop_z_below 0")
  las_check(las)
  st_crs(las) <- 25833
  
  # Add 50 cm buffer around the buildings
  buildings <- st_buffer(buildings, 0.5)
  
  # Normalization (subtract the DTM)
  gnd <- filter_ground(las)
  dtm <- rasterize_terrain(las, 1, knnidw())
  nlas <- normalize_height(las, knnidw())
  
  # Get the inverse of the buildings and crop it
  las_bbox <- st_as_sfc(st_bbox(las))
  buildings_aoi <- st_crop(buildings, las_bbox)
  erase_feature <- st_difference(las_bbox, st_union(buildings_aoi))

  las_clip <- clip_roi(nlas, erase_feature)
  
  return(las_clip)
}

# Folder paths
las_folder_path <- "C:/tree-canopy/data/400_5816/LAS"
buildings_path <- "C:/tree-canopy/data/400_5816/buildings.gpkg"
output_folder_path <- "C:/tree-canopy/data/400_5816/LAS_no_buildings"
dir.create(output_folder_path)

# Loop through all LAS files in the folder and remove buildings
las_files <- list.files(las_folder_path, pattern = "\\.las$", full.names = TRUE)
for (las_path in las_files) {
  cat("Processing:", las_path, "\n")
  las_nobuild <- remove_buildings_from_las(buildings_path, las_path)
  
  # Save the resulting LAS file without buildings
  output_file_name <- paste0(file_path_sans_ext(basename(las_path)), "_nobuild.las")
  output_file_path <- file.path(output_folder_path, output_file_name)
  writeLAS(las_nobuild, output_file_path)
}

