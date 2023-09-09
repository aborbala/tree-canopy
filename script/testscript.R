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

calculate_ratio_df <- function(df) {
  df %>%
    rowwise() %>%
    mutate(width_to_height_ratio = calculate_ratio(geometry))
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


#las_nobuild_path <- "../data/test 382_5826_ai/las_clipped.las"
las_nobuild_path <- "G:/Meine Ablage/data/382_5826_1/LAS_no_buildings/3dm_33_382_5827_1_be_nobuild.las"
#tif_directory <- "../data/test 382_5826_ai/382_5826_1_imgs/"
tif_directory <- "G:/Meine Ablage/data/382_5826_1/sliced_imgs_2020S"

las_nobuild <- readLAS(las_nobuild_path)
las_clip <- filter_poi(las_nobuild, Z >= 5)

files <- list.files(path = tif_directory, pattern = "\\.tif$", full.names = TRUE, recursive = FALSE)

#test <- files[75]
#print(test)
test <- "G:/Meine Ablage/data/382_5826_1/sliced_imgs_2020S/3dm_33_382_5827_1_be_0_2.tif"

ras <- rast(test)
#ras <- st_as_sf(read_stars(test))
st_crs(ras) <- 25833
ext <- st_bbox(ras)

##extract_crowns
las_aoi <- clip_roi(las_clip, ext)

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
plot(chm_p2r_05)

kernel <- matrix(1,3,3)
# Calculate focal ("moving window") values for each cell: smoothing steps with a median filter
chm_p2r_05_smoothed <- terra::focal(chm_p2r_05, w = kernel, fun = median, na.rm = TRUE)
ttops_chm_p2r_05_smoothed <- locate_trees(chm_p2r_05_smoothed, lmf((ws=10))) # using a greater window size
#par(mfrow=c(1,2))
plot(chm_p2r_05, col = height.colors(50))
plot(sf::st_geometry(ttops_chm_p2r_05_smoothed), add = TRUE, pch = 3)

### V1. Dalponte
algo <- dalponte2016(chm_p2r_05_smoothed, ttops_chm_p2r_05_smoothed)
las_tree <- segment_trees(las_aoi, algo)
crowns <- crown_metrics(las_tree, func = .stdtreemetrics, geom = "concave")
# Filter out small convex hulls, consider crowns over 1m2
crowns <- crowns[crowns$convhull_area > 10,]
plot(crowns["treeID"], main = "treeID")
#st_write(crowns, "crowns_dalponte.geojson")

crowns <- calculate_ratio_df(crowns)

crowns <- crowns[crowns$width_to_height_ratio > 0.5 & crowns$width_to_height_ratio < 1.5 ,]
plot(chm_p2r_05_smoothed, col = height.colors(50))
#plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
plot(crowns["convhull_area"], col = NA, border=2, add = TRUE)

### V2. Silva
terra::plot(rast(test))
#plot(chm_p2r_05_smoothed, col = height.colors(50), add = TRUE)

ttops <- locate_trees(chm_p2r_05_smoothed, lmf(10))
plot(ttops, add = TRUE)
algo_silva <- silva2016(chm_p2r_05_smoothed, ttops, ws = 10)
las_silva <- segment_trees(las_aoi, algo_silva)
ccm = ~custom_crown_metrics(z = Z, i = Intensity)
crowns_silva <- crown_metrics(las_silva, func = ccm, geom = "concave")
crowns_silva <- crowns_silva[crowns_silva$z_sd > 0.5,]
crowns_silva <- st_simplify(crowns_silva, dTolerance = 0.3)

#crowns_silva <- crowns_silva[crowns_silva$convhull_area > 10,]
#crowns_silva <- calculate_ratio_df(crowns_silva)
#crowns_silva <- crowns_silva[crowns_silva$width_to_height_ratio > 0.5 & crowns_silva$width_to_height_ratio < 1.5 ,]

#st_write(crowns_silva, "crowns_silva.geojson", delete_dsn = T)

plot(sf::st_geometry(ttops), add = TRUE, pch = 3, col=5)
plot(crowns_silva["treeID"], col = 'yellow', border=2, add=T)

###Watershed
algo_ws <- watershed(chm_p2r_05_smoothed)
las_ws <- segment_trees(las_aoi, algo_ws)
crowns_ws <- crown_metrics(las_ws, func = .stdtreemetrics, geom = "convex")
plot(crowns_ws["convhull_area"], col = pastel.colors(200))


