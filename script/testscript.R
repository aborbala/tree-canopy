# Load necessary libraries
library("rstudioapi")
library("Rcpp")
library("sf")
library("lidR")
library("stars")
library("terra")
library("raster")
library("dplyr")
library(viridisLite)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(patchwork)
library(showtext)
library(sysfonts)
library(terrainr)
library(EBImage)


font_add_google("Montserrat", "montserrat")
showtext_auto()

mako_colors <- viridis(n = 7, option = "mako")

# Define color variables
theme.main <- mako_colors[5]
theme.secondary <- mako_colors[2]
#color_red <- "#F21B1B"
#color_blue <- "#3647EB"
#color_green <- "#00C87F"
#color_lightblue <- "#00E3EC"
#color_yellow <- "#FFF200"
#color_magenta <- "magenta"
#color_orange <- "#fc9a19"
#color_purple <- "#9146fa"
#color_cyan <- "#4baabf"

color_red = "#ac3a38"        # Complementary red
color_blue = "#3870ac"       # Analogous blue
color_green = "#38ac74"      # Analogous green
color_lightblue = "#38aaac"  # Base color
color_yellow = "#FFF200"     # Triadic yellow-green
color_magenta = "#ac38aa"    # Triadic magenta
color_orange = "#fc9a19"     # Existing high-visibility orange
color_purple = "#9146fa"     # Existing high-visibility purple
color_cyan = "#4baabf"       # Existing high-visibility cyan

getwd()

## Helper Functions
# Calculate width-to-height ratio of a polygon
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

# Add width-to-height ratio to a data frame
calculate_ratio_df <- function(df) {
  df %>%
    rowwise() %>%
    mutate(width_to_height_ratio = calculate_ratio(geometry))
}

# Convert spatial data to dataframe format (ensuring EPSG:25833)
convert_sf_to_df <- function(sf_obj) {
  sf_obj <- st_transform(sf_obj, crs = 25833)  # Ensure correct CRS
  df <- as.data.frame(st_coordinates(sf_obj))  # Extract x, y
  df <- df[complete.cases(df), ]  # Remove NAs
  return(df)
}

# Convert polygon to dataframe
convert_polygon_to_df <- function(sf_obj) {
  sf_obj <- st_transform(sf_obj, crs = 25833)  # Ensure correct CRS
  poly_df <- st_coordinates(sf_obj) %>%
    as.data.frame() %>%
    rename(x = X, y = Y) %>%
    mutate(polygon_id = interaction(L1, L2, drop = TRUE))  # Ensure separate polygons
  return(poly_df)
}

calculate_metrics <- function(detected_trees, method_name) {
  intersections <- st_intersects(detected_trees, tree_cadaster_buffer, sparse = FALSE)
  true_positives <- sum(rowSums(intersections) > 0)  # Trees intersecting with cadaster
  false_positives <- nrow(detected_trees) - true_positives  # Trees outside cadaster buffer
  
  data.frame(
    Method = method_name,
    True_Positives = true_positives,
    False_Positives = false_positives
  )
}

# Custom crown metrics function
custom_crown_metrics <- function(z, i) { 
  metrics <- list(
    # Height-Based Metrics
    z_max = as.numeric(max(z, na.rm = TRUE)),        # Maximum height
    z_mean = as.numeric(mean(z, na.rm = TRUE)),      # Mean height
    z_sd = as.numeric(sd(z, na.rm = TRUE)),          # Height standard deviation
    z_range = as.numeric(max(z, na.rm = TRUE) - min(z, na.rm = TRUE)),  # Height range
    
    # LiDAR Intensity Metrics
    i_mean = as.numeric(mean(i, na.rm = TRUE)),      
    i_max = as.numeric(max(i, na.rm = TRUE)),        
    i_sd = as.numeric(sd(i, na.rm = TRUE)),          
    i_median = as.numeric(median(i, na.rm = TRUE)),  
    i_ratio = as.numeric(ifelse(max(i, na.rm = TRUE) > 0, mean(i, na.rm = TRUE) / max(i, na.rm = TRUE), NA))  
  )
  return(metrics)
}


# Main
## Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dop_tile <- "386_5818"

## Set up the theme
theme_set(
  theme_ipsum_rc(base_family = "montserrat") + 
    theme(
      panel.grid = element_blank(),      
      panel.background = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(margin = margin(t = 5), vjust = 1),
      axis.text.y = element_text(margin = margin(r = 5), hjust = 1),
      panel.grid.minor = element_blank()
    )
)

## Dynamically construct the paths
las_nobuild_path <- sprintf("G:/Meine Ablage/data/%s/LAS_no_buildings/3dm_33_386_5819_1_be_nobuild.las", dop_tile)
tif_directory <- sprintf("G:/Meine Ablage/data/%s/sliced_imgs_2020S", dop_tile)
test <- sprintf("G:/Meine Ablage/data/%s/sliced_imgs_2020S/3dm_33_386_5819_1_be_nobuild_8_5.tif", dop_tile)
tree_cadaster_path <- "C:/tree-canopy/data/386_5818/GRIS/tree_cadaster_3dm_33_386_5819_1_be_nobuild_8_5.gpkg"
tree_cadaster <- st_read(tree_cadaster_path)
st_crs(tree_cadaster) <- 25833

## Load and Filter LiDAR Data
las_nobuild <- readLAS(las_nobuild_path)

## Load Raster Data
ras <- rast(test)
ras <- st_as_sf(read_stars(test))
st_crs(ras) <- 25833
ext <- st_bbox(ras)

# Clip LiDAR Data to Raster Extent
las_unfiltered <- clip_roi(las_nobuild, ext)
las_above5 <- filter_poi(las_unfiltered, Z >= 5)

if (dim(las_above5@data)[1] < 5){
  return()
}
# -----------------------------------------------
### Canopy Height Model (CHM)
# -----------------------------------------------

## Point to raster
# algorithm for digital surface model computation based on a points-to-raster method: 
# for each pixel of the output raster the function attributes the height 
# of the highest point found. The subcircle tweak replaces each point 
# with 8 points around the original one.
# Points-to-raster algorithm with a resolution of 0.5 meters replacing each
# point by a 20 cm radius circle of 8 points
#chm_p2r_05 <- rasterize_canopy(las_unfiltered, 0.5, p2r(subcircle = 0.2, na.fill = tin()), pkg = "terra")
chm_p2r_05 <- rasterize_canopy(las_unfiltered, 0.5, p2r(subcircle = 0, na.fill = tin()), pkg = "terra")

## Triangulation
chm_dsmtin_base <- rasterize_canopy(las_unfiltered, res = 0.5, algorithm = dsmtin())
las_above3 <- filter_poi(las_unfiltered, Z >= 3)
chm_dsmtin_3 <- rasterize_canopy(las_above3, res = 0.5, algorithm = dsmtin(max_edge = 1.5))
chm_dsmtin_5 <- rasterize_canopy(las_above5, res = 0.5, algorithm = dsmtin(max_edge = 1.5))
las_above10 <- filter_poi(las_above5, Z >= 10)
chm_dsmtin_10 <- rasterize_canopy(las_above10, res = 0.5, algorithm = dsmtin(max_edge = 1.5))

# Pit-free
chm_pitfree <- rasterize_canopy(las_unfiltered, res = 0.5, pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 0.5)))
chm_pitfree_subcirlce <- rasterize_canopy(las_unfiltered, res = 0.5, pitfree( thresholds = c(0, 2, 5, 10, 15), subcircle = 0.15))

# Standard plotting
layers <- c(chm_p2r_05, chm_dsmtin_base, chm_dsmtin_5, chm_pitfree_subcirlce)
names(layers) <- c("chm_p2r_05", "chm_dsmtin_base", "chm_dsmtin_5", "chm_pitfree_subcirlce")
plot(layers)

# GGPLOT CHM methods
# Convert each raster to a dataframe
chm_p2r_05_df <- as.data.frame(chm_p2r_05, xy = TRUE)
chm_dsmtin_base_df <- as.data.frame(chm_dsmtin_base, xy = TRUE)
chm_dsmtin_5_df <- as.data.frame(chm_dsmtin_5, xy = TRUE)
chm_dsmtin_10_df <- as.data.frame(chm_dsmtin_10, xy = TRUE)
chm_pitfree_default_df <- as.data.frame(chm_pitfree, xy = TRUE)
chm_pitfree_df <- as.data.frame(chm_pitfree_subcirlce, xy = TRUE)

create_chm_plot <- function(df, title) {
  ggplot(df, aes(x = x, y = y, fill = Z)) +
    geom_tile() +
    scale_fill_viridis_c(option = "mako") +
    labs(title = title, fill = "Z", x = "Easting (m)", y = "Northing (m)") +
    scale_x_continuous(breaks = function(x) c(min(x), max(x))) +
    scale_y_continuous(breaks = function(y) c(min(y), max(y))) +
    coord_equal()
}

# Generate plots
p1 <- create_chm_plot(chm_p2r_05_df, "Point to raster")
p2 <- create_chm_plot(chm_dsmtin_base_df, "Triangulation")
p3 <- create_chm_plot(chm_dsmtin_5_df, "Triangulation over 5 m")
p4 <- create_chm_plot(chm_dsmtin_10_df, "Triangulation")
p5 <- create_chm_plot(chm_pitfree_default_df, "Pit-free")
p6 <- create_chm_plot(chm_pitfree_df, "Pit-free subcircle")

# Arrange plots in patchwork
chm_plot <- wrap_plots(p1, p2, p5, p6, ncol = 2) + 
  plot_layout(guides = 'collect') &  # Aligns legends and reduces spacing
  theme(plot.margin = margin(0, 0, 0, 0))  # Removes large gaps

# Print the combined plot
print(chm_plot)

# 🚀 Save with larger fonts
#ggsave("plots/chm_plot.png", plot = chm_plot, width = 12, height = 8, dpi = 360, scale = 2)
# -----------------------------------------------
### Locate trees
# -----------------------------------------------

# Local Maximum Filter with variable windows size
# Function for Deciduous Trees
ws_deciduous <- function(H) {
  return(3.09632 + 0.00895* (H^2))
}

# Function for Combined (Mixed) Trees
ws_combined <- function(H) {
  return(2.51503 + 0.00901 * (H^2))
}

heights <- seq(-5,30,0.5)
ws <- ws_deciduous(heights)
plot(heights, ws, type = "l",  ylim = c(0,10))

# Basemap
plot(chm_pitfree_subcirlce, col = magma(50))

## Cadaster
plot(sf::st_geometry(tree_cadaster),  add = TRUE, pch = 18, col="#F21B1B", cex = 3) #red

## LIDAR
ttops_LAS_ws_8 <- locate_trees(las_unfiltered, lmf(ws = 8, hmin=5))
plot(sf::st_geometry(ttops_LAS_ws_8), add = TRUE, pch = 1, col=color_blue, cex = 2, lwd = 3) #blue
ttops_LAS_ws_deciduous <- locate_trees(las_unfiltered, lmf(ws_deciduous, hmin=5, shape = c("square"), ws_args = list("Z")))
plot(sf::st_geometry(ttops_LAS_ws_deciduous), add = TRUE, pch = 2, col=color_green, cex = 2, lwd = 3) #blue
ttops_LAS_ws_combined <- locate_trees(las_unfiltered, lmf(ws_combined, hmin=5,  ws_args = list("Z")))
plot(sf::st_geometry(ttops_LAS_ws_combined), add = TRUE, pch = 3, col=color_red, cex = 2, lwd = 3) #blue

## CHM
# TIN
ttops_dsmtin_base <- locate_trees(chm_dsmtin_base, lmf(ws = 8, hmin=5))
plot(sf::st_geometry(ttops_dsmtin_base), add = TRUE, pch = 1, col=color_green, cex = 2, lwd = 3) #lightblue
ttops_dsmtin_deciduous <- locate_trees(chm_dsmtin_base, lmf(ws_deciduous, hmin=5))
plot(sf::st_geometry(ttops_dsmtin_deciduous), add = TRUE, pch = 2, col=color_yellow, cex = 2, lwd = 3) #green

# pit-free
ttops_pitfree_subcirlce_combined  <- locate_trees(chm_pitfree_subcirlce, lmf(ws = ws_combined, hmin = 5,  ws_args = list("Z")))
plot(sf::st_geometry(ttops_pitfree_subcirlce_combined), add = TRUE, pch = 2, col=color_yellow, cex = 2, lwd = 3) #green
ttops_pitfree_subcirlce_deciduous <- locate_trees(chm_pitfree_subcirlce, lmf(ws = ws_deciduous, hmin = 5,  ws_args = list("Z")))
plot(sf::st_geometry(ttops_pitfree_subcirlce_deciduous), add = TRUE, pch = 3, col=color_red, cex = 2, lwd = 3) #green

## compared to the Cadaster, these 2 extra tree tops are not even trees
# CHM is faster because there are less data to process, 
# but it is also more complex because the output depends on how the CHM has been built.

## GGPLOT Tree location
chm_df <- as.data.frame(chm_pitfree_subcirlce, xy = TRUE)
names(chm_df)[3] <- "value"

# Get first and last coordinate values for axis breaks
x_breaks <- range(chm_df$x, na.rm = TRUE)
y_breaks <- range(chm_df$y, na.rm = TRUE)

tree_cadaster_df <- convert_sf_to_df(tree_cadaster)
ttops_LAS_df <- convert_sf_to_df(ttops_LAS_ws_deciduous)
ttops_dsmtin_df <- convert_sf_to_df(ttops_dsmtin_deciduous)
ttops_pitfree_subcirlce_df <- convert_sf_to_df(ttops_pitfree_subcirlce_deciduous)
ttops_pitfree_subcirlce_combined_df <- convert_sf_to_df(ttops_pitfree_subcirlce_combined)
                                                      
# Base ggplot layer (CHM + Cadaster Trees) in UTM EPSG:25833
base_map <- ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(option = "mako", name = "CHM") +
  geom_point(data = tree_cadaster_df, aes(x = X, y = Y), color = color_orange, size = 4, shape = 18) +  # Convert to points
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_equal() +
  labs(x = "Easting (m)", y = "Northing (m)")

# Compare cadaster vs LIDAR detected trees
plot_lidar <- base_map +
  geom_point(data = ttops_LAS_df, aes(x = X, y = Y), color = color_yellow, size = 1, shape = 3,  stroke = 2) +
  labs(title = "Cadaster vs LiDAR Trees")

# Compare cadaster vs CHM detected trees (pitfree), ws_combined
plot_chm_3 <- base_map +
  geom_point(data = ttops_pitfree_subcirlce_combined_df, aes(x = X, y = Y), color = color_yellow, size = 1, shape = 3,  stroke = 2) +
  labs(title = "Cadaster vs CHM Trees, ws = combined")

# Compare cadaster vs CHM detected trees (pitfree), ws_decidous
plot_chm_5 <- base_map +
  geom_point(data = ttops_pitfree_subcirlce_df, aes(x = X, y = Y, color = "Detected Tree", shape = "Detected Tree"), size = 1, stroke = 2) +
  geom_point(data = tree_cadaster_df, aes(x = X, y = Y, color = "Cadaster Tree", shape = "Cadaster Tree"), size = 4) +
  scale_fill_viridis_c(option = "mako", name = "Height (m)") + 
  scale_color_manual(values = c("Cadaster Tree" = color_orange, "Detected Tree" = color_yellow), name = "Tree Type") +
  scale_shape_manual(values = c("Cadaster Tree" = 18, "Detected Tree" = 3), name = "Tree Type") +
  labs(title = "Cadaster vs CHM Trees, ws = deciduous") +
  theme(legend.position = "right") 

# Arrange all plots with patchwork (3 columns layout)
treetop_plot <- (plot_lidar | plot_chm_3 | plot_chm_5)

# Print the fixed combined plot
print(treetop_plot)
#ggsave("plots/treetop_plot.png", plot = treetop_plot, width = 12, height = 8, dpi = 360, scale = 1)

# ---------------------------------------------

## Postprocessing CHM: try to remove traffic signs and lamps
chm_pitfree_subcirlce[chm_pitfree_subcirlce < 5] <- NA
# Create a binary raster (e.g., threshold > 0 for tree areas)
binary_chm <- chm_pitfree_subcirlce > 5

# Identify connected components (clusters)
patches_chm <- patches(binary_chm, directions = 8, zeroAsNA = FALSE)

# Calculate patch sizes
patch_sizes <- freq(patches_chm) # Frequency table of patch IDs
small_patches <- patch_sizes$value[patch_sizes$count <= 5] # IDs of small patches

# Create a mask for small patches
small_patches_mask <- patches_chm %in% small_patches

# Remove small patches by masking them out
chm_pitfree_subcirlce_cleaned <- mask(chm_pitfree_subcirlce, small_patches_mask, maskvalue = TRUE)

# Set small patches to a unique value for visualization
#highlighted_patches <- ifel(small_patches_mask, 1, NA) # 1 for small patches, NA otherwise

# Plot the original raster
#plot(chm_dsmtin_5, main = "Original CHM with Highlighted Small Patches")

# Overlay small patches in magenta
#plot(highlighted_patches, col = "magenta", legend = FALSE, add = TRUE)

ttops_pitfree_subcirlce_cleaned <- locate_trees(chm_pitfree_subcirlce_cleaned, lmf(ws_deciduous, hmin=5, shape = c("square")))

terra::plot(rast(test))
plot(sf::st_geometry(ttops_pitfree_subcirlce_cleaned), add = TRUE, pch = 3,lwd = 3, col=color_yellow) #yellow

## GGPLOT Treetops on satellite img
# Convert raster to dataframe for ggplot
raster_raw_df <- as.data.frame(rast(test), xy = TRUE)
names(raster_raw_df)[3] <- "value"  # Standardize raster value column name

# Convert `ttops_dsmtin_5_cleaned` to dataframe
ttops_pitfree_subcirlce_cleaned_df <- as.data.frame(st_coordinates(ttops_pitfree_subcirlce_cleaned))  # Extract x/y

# Get first and last coordinate values for axis breaks
x_breaks <- range(raster_raw_df$x, na.rm = TRUE)
y_breaks <- range(raster_raw_df$y, na.rm = TRUE)

# Create ggplot with Raw Raster + Points
plot_raster_raw <- ggplot() +
  geom_spatial_rgb(
    data = test,
    mapping = aes(
      x = x,
      y = y,
      r = red,
      g = green,
      b = blue
    )
  ) +
  scale_fill_identity(guide = "none") +  # ✅ Use the exact raster colors without modification
  geom_point(data = ttops_pitfree_subcirlce_cleaned_df, aes(x = X, y = Y), 
             color = color_yellow, size = 2, shape = 3, stroke = 2) +  # ✅ Yellow cross points
  scale_y_continuous(breaks = y_breaks) +
  labs(title = "Tree tops (postprocessed)", x = "Easting (m)", y = "Northing (m)") +
  coord_equal() 

print(plot_raster_raw)
#ggsave("plots/treetop_postprocessed_satellite_plot.png", plot = plot_raster_raw, width = 12, height = 8, dpi = 100, scale = 1)

# ----------------------------------------------

## Get the stats of treetop detection
tree_cadaster <- tree_cadaster %>%
  mutate(buffer_size = ws_deciduous(baumhoehe))

# Create buffer zones
tree_cadaster_buffer <- st_buffer(tree_cadaster, dist = tree_cadaster$buffer_size)

plot(chm_pitfree_subcirlce_cleaned, col = mako(50))
plot(st_geometry(tree_cadaster_buffer), col = "transparent", border = color_red, add = T, main = "Tree Buffers Based on Height")
plot(st_geometry(tree_cadaster), col = color_red, pch = 16, add = TRUE)
# GGPLOT
# Convert raster to data frame for ggplot
#chm_df <- as.data.frame(chm_pitfree_subcirlce_cleaned, xy = TRUE)
#x_breaks <- range(chm_df$x, na.rm = TRUE)
#y_breaks <- range(chm_df$y, na.rm = TRUE)
#names(chm_df)[3] <- "value"

# Convert tree buffer and tree cadaster to sf objects if not already
tree_cadaster_buffer_df <- convert_sf_to_df(tree_cadaster_buffer)
tree_cadaster_buffer_sf <- st_as_sf(tree_cadaster_buffer)  # Preserve buffer as sf POLYGON
tree_cadaster_buffer_poly <- st_coordinates(tree_cadaster_buffer_sf) %>%
  as.data.frame() %>%
  rename(x = X, y = Y) %>%
  mutate(polygon_id = interaction(L1, L2, drop = TRUE))  # Correctly separates polygons and their rings


# Create the base map
base_map <- ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = tree_cadaster_buffer_poly, 
               aes(x = x, y = y, group = polygon_id), 
               color = "red", fill = NA, linewidth = 0.7) +
  geom_point(data = tree_cadaster_df, aes(x = X, y = Y), color = color_red, size = 2) +
  scale_fill_viridis_c(option = "mako", name = "CHM") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  labs(x = "Easting (m)", y = "Northing (m)", title = "Tree Buffers Based on Height")

# Display the plot
print(base_map)
## -------------------------------
calculate_iou_per_tree <- function(detected_trees) {
  detected_trees <- detected_trees %>%
    mutate(buffer_size = ws_deciduous(Z))  # Assign buffer based on height
  
  detected_tree_buffer <- st_buffer(detected_trees, dist = detected_trees$buffer_size)
  
  iou_values <- sapply(1:nrow(detected_trees), function(i) {
    detected_tree <- detected_tree_buffer[i, ]  # Extract single detected tree
    nearest_cadaster <- st_nearest_feature(detected_tree, tree_cadaster_buffer)  # Find closest true tree
    matched_cadaster <- tree_cadaster_buffer[nearest_cadaster, ]
    
    intersections <- st_intersection(detected_tree, matched_cadaster)
    if (nrow(intersections) == 0) {
      return(0)  # No intersection = IoU of 0
    } else {
      intersection_area <- sum(as.numeric(st_area(intersections)))
      union_area <- sum(as.numeric(st_area(st_union(detected_tree, matched_cadaster))))
      return(ifelse(union_area > 0, intersection_area / union_area, 0))
    }
  })
}
# Function to calculate true positives and false positives
calculate_metrics <- function(detected_trees, method_name) {
  # Create buffer around detected trees based on height
  detected_trees <- detected_trees %>%
    mutate(buffer_size = ws_deciduous(Z))
  detected_tree_buffer <- st_buffer(detected_trees, dist = detected_trees$buffer_size)
  
  # PER-TREE: Compute intersection and union areas (convert units to numeric)
  mean_IoU <- mean(calculate_iou_per_tree(detected_trees), na.rm = TRUE)
  median_IoU <- median(calculate_iou_per_tree(detected_trees), na.rm = TRUE)

  # GLOBAL: Compute Global IoU by dissolving geometries
  detected_tree_buffer_dissolved <- st_union(detected_tree_buffer)
  tree_cadaster_buffer_dissolved <- st_union(tree_cadaster_buffer)
  
  global_intersection_area <- sum(as.numeric(st_area(st_intersection(detected_tree_buffer_dissolved, tree_cadaster_buffer_dissolved))))
  global_union_area <- sum(as.numeric(st_area(st_union(detected_tree_buffer_dissolved, tree_cadaster_buffer_dissolved))))
  global_IoU <- ifelse(global_union_area > 0, global_intersection_area / global_union_area, 0)
  
  # True Positives (detections overlapping with cadaster buffer)
  intersections <- st_intersects(detected_trees, tree_cadaster_buffer, sparse = FALSE)
  true_positives <- sum(rowSums(intersections) > 0)
  
  # False Positives (detected trees not intersecting with cadaster)
  false_positives <- nrow(detected_trees) - true_positives
  # False Negatives (cadaster trees not detected)
  false_negatives <- nrow(tree_cadaster) - true_positives
  
  # Precision, Recall, and F1 Score
  precision <- ifelse((true_positives + false_positives) > 0, true_positives / (true_positives + false_positives), 0)
  recall <- ifelse((true_positives + false_negatives) > 0, true_positives / (true_positives + false_negatives), 0)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
  
  # Deviation (Absolute & Relative)
  count_cadaster <- nrow(tree_cadaster)
  count_detected <- nrow(detected_trees)
  absolute_deviation <- count_detected - count_cadaster
  relative_deviation <- ifelse(count_cadaster > 0, (count_detected / count_cadaster) * 100, NA)
  

  data.frame(
    Method = method_name,
    True_Positives = true_positives,
    False_Positives = false_positives,
    False_Negatives = false_negatives,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score,
    mean_IoU = mean_IoU,
    median_IoU = median_IoU,
    Global_IoU = global_IoU,
    Count_Cadaster = count_cadaster,
    Count_Detected = count_detected,
    Absolute_Deviation = absolute_deviation,
    Relative_Deviation = relative_deviation
    
  )
}

# Calculate metrics for each method
results <- bind_rows(
  calculate_metrics(ttops_LAS_ws_8, "LIDAR ws=8"),
  calculate_metrics(ttops_LAS_ws_deciduous, "LIDAR ws=d"),
  calculate_metrics(ttops_dsmtin_deciduous, "CHM dsmtin ws=d"),
  calculate_metrics(ttops_pitfree_subcirlce_deciduous, "CHM pitfree ws=d"),
  calculate_metrics(ttops_pitfree_subcirlce_cleaned, "CHM pitfree ws=d [cleaned]")
)


# Sort results by True Positives in descending order
#results_sorted <- results %>%
#  arrange(desc(True_Positives))

# Convert Method to a factor with levels in sorted order
results$Method <- factor(results$Method, levels = results$Method)

# Convert data to long format
results_long <- results %>%
  pivot_longer(cols = c(mean_IoU, Global_IoU), 
               names_to = "IoU_Type", 
               values_to = "IoU_Score") %>%
  mutate(IoU_Type = recode(IoU_Type, "mean_IoU" = "Per-Tree IoU", "Global_IoU" = "Global IoU"))

# IoU bar per tree vs Global
iou_plot <- ggplot(results_long, aes(x = Method, y = IoU_Score, fill = IoU_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Per-Tree IoU" = theme.secondary, "Global IoU" = theme.main)) +  # Improved colors
  labs(title = "Comparison of IoU Metrics",
       x = "Method",
       y = "IoU Score",
       fill = "IoU Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),  # Subtle grid
    panel.grid.minor = element_blank()
  )

# Precesion-Recall plot
precision_recall_ttops_plot <- ggplot(results, aes(x = Recall, y = Precision, label = Method, color = F1_Score)) +
  geom_point(size = 5) + 
  geom_text(color = "black", size = 4, vjust = -1, hjust = 0.5, family = "montserrat") +  # Ensuring font consistency
  scale_color_viridis_c(option = "mako", begin = 0.2, end = 0.9, name = "F1 Score") +  # Use "mako" with adjusted range
  labs(title = "Precision-Recall Tradeoff",
       x = "Recall",
       y = "Precision") +
  theme(
    legend.position = "right",
    legend.margin = margin(l = 30)
  ) +
  coord_equal(ratio=2,clip = "off")

# Arrange all plots side by side
eval_plot <-  iou_plot | precision_recall_ttops_plot

# Display the final combined plot
print(eval_plot)
print(precision_recall_ttops_plot)
#ggsave("plots/precision_recall_ttops_plot.png", plot = precision_recall_ttops_plot, width = 12, height = 8, dpi = 100, scale = 1)

# -------------------------------------
# Count of trees detected
# Calculate total detected tree count for each method
total_counts <- data.frame(
  Method = c("LIDAR ws=8", "CHM pitfree cleaned ws=d", "CHM pitfree ws=d","CHM dsmtin ws=d", "LIDAR ws=d"),
  Detected_Trees = c(nrow(ttops_LAS_ws_8), 
                     nrow(ttops_pitfree_subcirlce_cleaned), 
                     nrow(ttops_pitfree_subcirlce_deciduous),
                     nrow(ttops_dsmtin_deciduous),
                     nrow(ttops_LAS_ws_deciduous))
)

# Get the cadaster tree count
cadaster_count <- nrow(tree_cadaster)

# Compute absolute and percentage deviation
total_counts <- total_counts %>%
  mutate(
    Absolute_Deviation = Detected_Trees - cadaster_count,
    Percentage_Deviation = (Absolute_Deviation / cadaster_count) * 100
  ) %>%
  arrange(desc(Absolute_Deviation))  # Sort in descending order

# Convert Method to a factor with sorted levels
total_counts$Method <- factor(total_counts$Method, levels = total_counts$Method)

# Plot: ABSOLUTE DEVIATION: Detected Trees vs. Cadaster Baseline (Sorted)
ggplot(total_counts, aes(x = Method, y = Absolute_Deviation, fill = Absolute_Deviation > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Tree Count: Deviation from Cadaster",
       x = "Detection Method",
       y = "Deviation from Cadaster",
       fill = "Deviation") +
  scale_fill_manual(values = c("TRUE" = theme.main, "FALSE" = theme.secondary)) +  # Green for overestimation, Red for underestimation
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot: Percentage Deviation
plot_percentage_deviation <- ggplot(total_counts, aes(x = Method, y = Percentage_Deviation, fill = Percentage_Deviation > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Tree Count: Percentage Deviation from Cadaster",
       x = "Detection Method",
       y = "Deviation (%)",
       fill = "Deviation") +
  scale_fill_manual(values = c("TRUE" = theme.main, "FALSE" = theme.secondary)) +  # Green for overestimation, Red for underestimation
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_text(angle = 90, hjust = 1))

# Arrange all plots side by side
eval_plot <-  plot_percentage_deviation

# Display the final combined plot
print(eval_plot)
#--------------------------------------

## OUTCOME:
# best lmf:
chm <- chm_pitfree_subcirlce_cleaned
# best chm vs lidar
ttops <- ttops_pitfree_subcirlce_cleaned

# -----------------------------------------------
### Tree Segmentation
# -----------------------------------------------
plot(chm)
plot(sf::st_geometry(ttops), add = TRUE, pch = 3, col=color_red, cex = 2, lwd = 3) #green
ccm = ~custom_crown_metrics(z = Z, i = Intensity)

## Dalponte Method (raster-based)
algo_dalponte <- dalponte2016(chm, ttops,   th_tree = 2,    # Minimum tree height
                              th_seed = 0.45, # Seed threshold for initial growth
                              th_cr = 0.65,   # Crown merging threshold
                              max_cr = 20    # Max crown diameter (20 pixels = 10m for 0.5m CHM)
)
las_dalponte <- segment_trees(las_unfiltered, algo_dalponte)

#tree10 <- filter_poi(las_tree, treeID == 11)
#plot(tree10, size = 8, bg = "white")
crowns_dalponte <- crown_metrics(las_dalponte, func = ccm, geom = "concave")
## Post processing: Filter out small convex hulls, consider crowns over 1m2
plot(crowns_dalponte["i_ratio"], main = "z_max")
#crowns <- calculate_ratio_df(crowns)
#crowns <- crowns[crowns$width_to_height_ratio > 0.5 & crowns$width_to_height_ratio < 1.5 ,]
plot(chm, col = mako(50))
#plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
plot(crowns_dalponte["z_max"], col = "transparent", border = color_cyan, lwd = 3, add = TRUE)
plot(sf::st_geometry(ttops), add = TRUE, pch = 3, col=color_red, cex = 2, lwd = 3)

## Silva method
#terra::plot(rast(test))
#ttops <- locate_trees(chm_p2r_05_smoothed, lmf(10))
algo_silva <- silva2016(chm, ttops, max_cr = 0.6, exclusion = 0.3)
las_silva <- segment_trees(las_unfiltered, algo_silva)
crowns_silva <- crown_metrics(las_silva, func = ccm, geom = "concave")
#crowns_silva <- crowns_silva[crowns_silva$z_sd > 0.5,]
#crowns_silva <- st_simplify(crowns_silva, dTolerance = 0.3)

#crowns_silva <- crowns_silva[crowns_silva$convhull_area > 10,]
#crowns_silva <- calculate_ratio_df(crowns_silva)
#crowns_silva <- crowns_silva[crowns_silva$width_to_height_ratio > 0.5 & crowns_silva$width_to_height_ratio < 1.5 ,]

#st_write(crowns_silva, "crowns_silva.geojson", delete_dsn = T)
plot(chm, col = magma(50))
#plot(crowns_silva["convhull_area"], col = NA, border=6, add = TRUE)
plot(crowns_silva["i_mean"], col="transparent", border = color_cyan, lwd = 3, add = TRUE)
plot(sf::st_geometry(ttops), add = TRUE, pch = 3, col=color_red,  cex = 2,lwd = 2)

## GGPLOT Segmentation
convert_polygon_to_df <- function(sf_obj) {
  sf_obj <- st_transform(sf_obj, crs = 25833)  # Ensure correct CRS
  poly_df <- st_coordinates(sf_obj) %>%
    as.data.frame() %>%
    rename(x = X, y = Y) %>%
    mutate(polygon_id = interaction(L1, L2, drop = TRUE))  # Ensure separate polygons
  return(poly_df)
}

crowns_dalponte_df <- convert_polygon_to_df(crowns_dalponte)
crowns_silva_df <- convert_polygon_to_df(crowns_silva)

# Base ggplot layer (CHM + Cadaster Trees) in UTM EPSG:25833
base_map <- ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(option = "mako", name = "Height (m)") +
  #geom_point(data = tree_cadaster_df, aes(x = X, y = Y), color = color_orange, size = 4, shape = 18) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_equal() +
  labs(x = "Easting (m)", y = "Northing (m)")


# Dalponte method visualization
# Base ggplot layer (CHM + Cadaster Trees) in UTM EPSG:25833
base_map <- ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(option = "mako", name = "Height (m)") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_equal() +
  labs(x = "Easting (m)", y = "Northing (m)")

# Dalponte method visualization
plot_dalponte <- base_map +
  geom_polygon(data = crowns_dalponte_df, aes(x = x, y = y, group = polygon_id),
               color = color_magenta, fill = NA, linewidth = 1) +
  geom_point(data = ttops_LAS_df, aes(x = X, y = Y), color = color_yellow, shape = 3, size = 1, stroke = 2) +
  geom_point(data = tree_cadaster_df, aes(x = X, y = Y), color = color_orange, size = 4, shape = 18) +
  labs(title = "Dalponte 2016 Segmentation")

# Silva method visualization (Fixing the legend issue)
plot_silva <- base_map +
  geom_polygon(data = crowns_silva_df, aes(x = x, y = y, group = polygon_id),
               color = color_magenta, fill = NA, linewidth = 1) +
  geom_point(data = tree_cadaster_df, aes(x = X, y = Y, color = "Cadaster Tree", shape = "Cadaster Tree"), size = 4, shape = 18) +
  geom_point(data = ttops_LAS_df, aes(x = X, y = Y, color = "Detected Tree", shape = "Detected Tree"),shape = 3, size = 1, stroke = 2) +
  scale_color_manual(name = "Tree Type", values = c("Cadaster Tree" = color_orange, "Detected Tree" = color_yellow)) +
  scale_shape_manual(name = "Tree Type", values = c("Cadaster Tree" = 18, "Detected Tree" = 3)) + 
  theme(legend.position = "right")+
  labs(title = "Silva 2016 Segmentation")


# Arrange both plots side by side
combined_plot <- plot_dalponte + plot_silva

# Print the fixed combined plot
print(combined_plot)
#ggsave("plots/dalponte_silva_plot.png", plot = combined_plot, width = 12, height = 8, dpi = 100, scale = 1)

# Evaluation
calculate_iou_for_crowns <- function(segmented_crowns, method_name, tree_cadaster_buffer) {
  # Compute intersection and union areas for per-crown IoU
  intersections <- st_intersection(segmented_crowns, tree_cadaster_buffer)
  intersection_area <- sum(as.numeric(st_area(intersections)), na.rm = TRUE)
  
  union_shapes <- st_union(st_union(segmented_crowns), st_union(tree_cadaster_buffer))
  union_area <- sum(as.numeric(st_area(union_shapes)), na.rm = TRUE)
  
  per_tree_iou <- ifelse(union_area > 0, intersection_area / union_area, 0)
  
  # Global IoU: dissolve all crowns and cadaster trees
  segmented_crowns_dissolved <- st_union(segmented_crowns)
  tree_cadaster_dissolved <- st_union(tree_cadaster_buffer)
  
  global_intersection_area <- sum(as.numeric(st_area(st_intersection(segmented_crowns_dissolved, tree_cadaster_dissolved))))
  global_union_area <- sum(as.numeric(st_area(st_union(segmented_crowns_dissolved, tree_cadaster_dissolved))))
  global_IoU <- ifelse(global_union_area > 0, global_intersection_area / global_union_area, 0)
  
  # Count overlap-based true positives
  matched_crowns <- st_intersects(segmented_crowns, tree_cadaster_buffer, sparse = FALSE)
  true_positives <- sum(rowSums(matched_crowns) > 0)
  
  # False positives: detected crowns that do not overlap with cadaster trees
  false_positives <- nrow(segmented_crowns) - true_positives
  
  # False negatives: cadaster trees not covered by any segmented crown
  false_negatives <- nrow(tree_cadaster_buffer) - true_positives
  
  # Compute precision, recall, and F1-score
  precision <- ifelse((true_positives + false_positives) > 0, true_positives / (true_positives + false_positives), 0)
  recall <- ifelse((true_positives + false_negatives) > 0, true_positives / (true_positives + false_negatives), 0)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
  
  # Count cadaster trees and detected crowns
  count_cadaster <- nrow(tree_cadaster_buffer)
  count_detected <- nrow(segmented_crowns)
  
  # Compute absolute and relative deviation
  absolute_deviation <- count_detected - count_cadaster
  relative_deviation <- ifelse(count_cadaster > 0, (count_detected / count_cadaster) * 100, NA)
  
  # Return results as a dataframe
  data.frame(
    Method = method_name,
    True_Positives = true_positives,
    False_Positives = false_positives,
    False_Negatives = false_negatives,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score,
    mean_IoU = per_tree_iou,
    Global_IoU = global_IoU,
    Count_Cadaster = count_cadaster,
    Count_Detected = count_detected,
    Absolute_Deviation = absolute_deviation,
    Relative_Deviation = relative_deviation
  )
}

results_dalponte <- calculate_iou_for_crowns(crowns_dalponte, "Dalponte", tree_cadaster_buffer)
results_silva <- calculate_iou_for_crowns(crowns_silva, "Silva", tree_cadaster_buffer)

# Combine results
crown_results <- rbind(results_dalponte, results_silva)

# Display results
print(crown_results)

# Convert data to long format
results_long <- crown_results %>%
  pivot_longer(cols = c(mean_IoU, Global_IoU), 
               names_to = "IoU_Type", 
               values_to = "IoU_Score") %>%
  mutate(IoU_Type = recode(IoU_Type, "mean_IoU" = "Per-Tree IoU", "Global_IoU" = "Global IoU"))


iou_plot <- ggplot(results_long, aes(x = Method, y = IoU_Score, fill = IoU_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Per-Tree IoU" = mako_colors[2], "Global IoU" = mako_colors[5])) +  # Improved colors
  labs(title = "Comparison of IoU Metrics",
       x = "Method",
       y = "IoU Score",
       fill = "IoU Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),  # Subtle grid
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 20, 10)
  )

# Precision-Recall plot (Removed coord_fixed and adjusted margins)
precision_recall_plot <- ggplot(crown_results, aes(x = Recall, y = Precision, label = Method, color = F1_Score)) +
  geom_point(size = 5) + 
  geom_text(color = "black", size = 4, vjust = -0.5, hjust = -0.1, family = "montserrat") +  # Ensuring font consistency
  scale_color_viridis_c(option = "mako", begin = 0.2, end = 0.9, name = "F1 Score") +  # Use "mako" with adjusted range
  labs(title = "Precision-Recall Tradeoff",
       x = "Recall",
       y = "Precision") +
  theme(
    legend.position = "right",
    legend.margin = margin(l = 30),
    plot.margin = margin(20, 10, 50, 10)  # Increased margin for better height
  )+
  coord_fixed(ratio = 5, clip= "off")

# Arrange all plots with different heights
segmentation_eval_plot <- iou_plot | precision_recall_plot +  
  plot_layout(heights = c(0.5, 5), widths = c(1, 1))  # Precision-Recall plot gets more vertical space

# Display the final combined plot
print(segmentation_eval_plot)
print(precision_recall_plot)
#ggsave("plots/dalponte_silva_precision_recall_plot.png", plot = precision_recall_plot, width = 12, height = 8, dpi = 100, scale = 1)

