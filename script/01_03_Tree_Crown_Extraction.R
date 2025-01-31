# Load necessary libraries
library("rstudioapi")
library("Rcpp")
library("sf")
library("lidR")
library("stars")
library("terra")
library("raster")
library("dplyr")

# Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Required file paths
# Extract crowns for small tif tiles
tif_directory <- "G:/Meine Ablage/data/386_5818/sliced_imgs_2020S/"

las_nobuild_path <- "C:/tree-canopy/data/386_5818/LAS_no_buildings"
las_files <- list.files(path = las_nobuild_path, pattern = "\\.las$", full.names = TRUE, recursive = FALSE)

# Outputs
crowns_path <- "G:/Meine Ablage/data/386_5818/crowns/"

### Calculate width to height ratio to eliminate elongated polygons
### Using bbox instead of oo-bbox
calculate_ratio <- function(polygon) {
  # Calculate the minimum bounding box
  mbr <- st_bbox(polygon)
  
  # The width and height of the MBR can be calculated from its coordinates
  width = mbr["xmax"] - mbr["xmin"]
  height = mbr["ymax"] - mbr["ymin"]
  
  # Calculate the width-to-height ratio
  ratio = width / height
  
  return(ratio)
}
custom_crown_metrics <- function(z, i) { # user-defined function
  metrics <- list(
    z_max = max(z),   # max height
    z_sd = sd(z),     # vertical variability of points
    i_mean = mean(i), # mean intensity
    i_max  = max(i)   # max intensity
  )
  return(metrics) # output
}

calculate_ratio_df <- function(df) {
  df %>%
    rowwise() %>%
    mutate(width_to_height_ratio = calculate_ratio(geometry))
}

### Function to extract crowns from Lidar data
extract_crowns <- function(las_clip, bbox) {
  print("extracting crowns...")
  las_aoi <- clip_roi(las_clip, bbox)

  # If no LAS point is present, skip
  if (dim(las_aoi@data)[1] < 10){
    print("extracting crowns... NO LAS present, return...")
    return()
  }
  print("extracting crowns... LAS present!")

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
  
  # Use dalponte algorithm
  #ttops_chm_p2r_05_smoothed <- locate_trees(chm_p2r_05_smoothed, lmf((ws=8))) # using a greater window size
  #algo <- dalponte2016(chm_p2r_05_smoothed, ttops_chm_p2r_05_smoothed)
  #las_tree <- segment_trees(las_aoi, algo)
  #crowns <- crown_metrics(las_tree, func = .stdtreemetrics, geom = "convex")
  #crowns <- crowns[crowns$convhull_area > 1,]
  
  # Use Silva algorithm
  ttops_chm_p2r_05_smoothed <- locate_trees(chm_p2r_05_smoothed, lmf(8))
  algo <- silva2016(chm_p2r_05_smoothed, ttops_chm_p2r_05_smoothed)
  las_tree <- segment_trees(las_aoi, algo)
  #crowns <- crown_metrics(las_tree, func = .stdtreemetrics, geom = "convex")
  #crowns <- crowns[crowns$convhull_area > 10,]
  #crowns <- calculate_ratio_df(crowns)
  #crowns <- crowns[crowns$width_to_height_ratio > 0.5 & crowns$width_to_height_ratio < 1.5 ,]
  ccm = ~custom_crown_metrics(z = Z, i = Intensity)
  crowns <- crown_metrics(las_tree, func = ccm, geom = "concave")
  crowns <- crowns[crowns$z_sd > 0.5,]
  #crowns <- calculate_ratio_df(crowns)
  #crowns <- crowns[crowns$width_to_height_ratio > 0.5 & crowns$width_to_height_ratio < 1.5 ,]
  crowns <- st_simplify(crowns, dTolerance = 0.3)
  
  return(crowns)
}

### Function to process all .tif files in a given directory
# process_tif_files <- function(dir_path, crowns_path) {
#   files <- list.files(path = dir_path, pattern = "\\.tif$", full.names = TRUE, recursive = FALSE)
#   
#   lapply(files, function(x) {
#     print(x)
#     ras <- st_as_sf(read_stars(x))
#     st_crs(ras) <- 25833
#     ext <- st_bbox(ras)
#     crowns <- extract_crowns(las_clip, ext)
#     
#     if (!is.null(crowns)) {
#       filename <- strsplit(sub(".*/", "", x), "\\.")[[1]][1]
#       st_write(crowns, paste(crowns_path, filename, ".geojson", sep = ""), delete_dsn = T)
#     }
#   })
# }
### Function to process all .tif files in a given directory
process_tif_files <- function(dir_path, crowns_path, las_files) {
  tif_files <- list.files(path = dir_path, pattern = "\\.tif$", full.names = TRUE, recursive = FALSE)
  lapply(tif_files, function(tif_file) {
    print(tif_file)
    ras <- st_as_sf(read_stars(tif_file))
    st_crs(ras) <- 25833
    ext <- st_bbox(ras)
    
    # Extract base name of the tif file
    tif_base_name <- strsplit(basename(tif_file), "_be_")[[1]][1]
    
    # Find the matching las file
    las_file <- las_files[grep(tif_base_name, las_files)]
    las_clip <- readLAS(las_file)
    
    # Only consider crowns taller than 5 m
    las_clip <- filter_poi(las_clip, Z >= 5)
    chm <- rasterize_canopy(las_clip, 1, p2r())
    
    crowns <- extract_crowns(las_clip, ext)
    
    if (!is.null(crowns)) {
      filename <- strsplit(sub(".*/", "", tif_file), "\\.")[[1]][1]
      st_write(crowns, paste(crowns_path, filename, ".geojson", sep = ""), delete_dsn = T)
    }
  })
}


### Processing

#las_nobuild <- readLAS(las_nobuild_path)

# Only consider crowns taller than 10 m
#las_clip <- filter_poi(las_nobuild, Z >= 5)
#chm <- rasterize_canopy(las_clip, 1, p2r())
#plot(chm)

#process_tif_files(tif_directory, crowns_path)
process_tif_files(tif_directory, crowns_path, las_files)

