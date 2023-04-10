# Load necessary libraries
library("rstudioapi")
library("Rcpp")
library("sf")
library("lidR")
library("stars")
library("terra")
library("raster")

# Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Required file paths
buildings_path <- "../data/test 382_5826_ai/buildings.gpkg"
las_path <- "../data/test 382_5826_ai/3dm_33_382_5826_1_be.las"
# Extract crowns for small tif tiles
tif_directory <- "../data/test 382_5826_ai/382_5826_1_imgs/"

# Outputs
las_nobuild_path <- "../data/test 382_5826_ai/las_clipped.las"
crowns_path <- "../data/test 382_5826_ai/382_5826_1_crowns/"

### Function to remove buildings from Lidar data
remove_buildings_from_las <- function(buildings_path, las_path) {
  buildings <- st_read(buildings_path)
  crs <- st_crs(buildings)
  las <- readLAS(las_path, filter = "-drop_z_below 0")
  las_check(las)
  st_crs(las) <- 25833
  
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

### Function to extract crowns from Lidar data
extract_crowns <- function(las_clip, bbox) {
  las_aoi <- clip_roi(las_clip, bbox)
  
  if (dim(las_aoi@data)[1] < 5){
    return()
  }
  
  # algorithm for digital surface model computation based on a points-to-raster method: 
  # for each pixel of the output raster the function attributes the height 
  # of the highest point found. The subcircle tweak replaces each point 
  # with 8 points around the original one.
  # Points-to-raster algorithm with a resolution of 0.5 meters replacing each
  # point by a 20 cm radius circle of 8 points
  chm_p2r_05 <- rasterize_canopy(las_aoi, 0.5, p2r(subcircle = 0.2), pkg = "terra")
  
  # Calculate focal ("moving window") values for each cell: smoothing steps with a median filter
  kernel <- matrix(1,3,3)
  chm_p2r_05_smoothed <- terra::focal(chm_p2r_05, w = kernel, fun = median, na.rm = TRUE)
  ttops_chm_p2r_05_smoothed <- locate_trees(chm_p2r_05_smoothed, lmf((ws=8))) # using a greater window size
  algo <- dalponte2016(chm_p2r_05_smoothed, ttops_chm_p2r_05_smoothed)
  las_tree <- segment_trees(las_aoi, algo)
  crowns <- crown_metrics(las_tree, func = .stdtreemetrics, geom = "convex")
  
  return(crowns)
}

### Function to process all .tif files in a given directory
process_tif_files <- function(dir_path, crowns_path) {
  files <- list.files(path = dir_path, pattern = "*.tif", full.names = TRUE, recursive = FALSE)
  
  lapply(files, function(x) {
    print(x)
    ras <- st_as_sf(read_stars(x))
    st_crs(ras) <- 25833
    ext <- st_bbox(ras)
    crowns <- extract_crowns(las_clip, ext)
    
    if (!is.null(crowns)) {
      filename <- strsplit(sub(".*/", "", x), "\\.")[[1]][1]
      st_write(crowns, paste(crowns_path, filename, ".geojson", sep = ""), delete_dsn = T)
    }
  })
}

### Processing

# Remove building points from LAS
las_nobuild <- remove_buildings_from_las(buildings_path, las_path)
writeLAS(las_nobuild, las_nobuild_path)

#las_nobuild <- readLAS(las_no_buildings_path)
# Only consider crowns taller than 10 m
las_clip <- filter_poi(las_nobuild, Z >= 10)
chm <- rasterize_canopy(las_clip, 1, pitfree(thr, edg))
plot(chm)

process_tif_files(tif_directory, crowns_path)
