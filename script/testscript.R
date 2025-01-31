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

## Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dop_tile <- "386_5818"
## Set up the theme
theme_set(
  theme_ipsum_rc(base_family = "montserrat") + 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),      
      panel.background = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(margin = margin(t = 5), vjust = 1),  # Move up
      axis.text.y = element_text(margin = margin(r = 5), hjust = 1)   # Move right
    )
)
theme.main <- "#009999"
theme.secondary <- "#6b2e34"

# Dynamically construct the paths
las_nobuild_path <- sprintf("G:/Meine Ablage/data/%s/LAS_no_buildings/3dm_33_386_5819_1_be_nobuild.las", dop_tile)
tif_directory <- sprintf("G:/Meine Ablage/data/%s/sliced_imgs_2020S", dop_tile)
test <- sprintf("G:/Meine Ablage/data/%s/sliced_imgs_2020S/3dm_33_386_5819_1_be_nobuild_8_5.tif", dop_tile)
tree_cadaster_path <- "C:/tree-canopy/data/386_5818/GRIS/tree_cadaster_3dm_33_386_5819_1_be_nobuild_8_5.gpkg"
tree_cadaster <- st_read(tree_cadaster_path)
st_crs(tree_cadaster) <- 25833

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
chm_p2r_05 <- rasterize_canopy(las_unfiltered, 0.5, p2r(subcircle = 0.2, na.fill = tin()), pkg = "terra")

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

# GGPLOT
# Convert each raster to a dataframe
chm_p2r_05_df <- as.data.frame(chm_p2r_05, xy = TRUE)
chm_dsmtin_base_df <- as.data.frame(chm_dsmtin_base, xy = TRUE)
chm_dsmtin_5_df <- as.data.frame(chm_dsmtin_5, xy = TRUE)
chm_dsmtin_10_df <- as.data.frame(chm_dsmtin_10, xy = TRUE)
chm_pitfree_df <- as.data.frame(chm_pitfree_subcirlce, xy = TRUE)

# Create individual plots for each raster
p1 <- ggplot(chm_p2r_05_df, aes(x = x, y = y, fill = Z)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Point to raster", fill = "Z") +
  scale_x_continuous(breaks = function(x) c(min(x), max(x))) +
  scale_y_continuous(breaks = function(y) c(min(y), max(y))) + 
  labs(x = "Easting (m)", y = "Northing (m)") +
  coord_equal()

p2 <- ggplot(chm_dsmtin_base_df, aes(x = x, y = y, fill = Z)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Triangulation over 0 m", fill = "Z") +
  scale_x_continuous(breaks = function(x) c(min(x), max(x))) +
  scale_y_continuous(breaks = function(y) c(min(y), max(y))) + 
  labs(x = "Easting (m)", y = "Northing (m)") +
  coord_equal()

p3 <- ggplot(chm_dsmtin_5_df, aes(x = x, y = y, fill = Z)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Triangulation over 5 m", fill = "Z") +
  scale_x_continuous(breaks = function(x) c(min(x), max(x))) +
  scale_y_continuous(breaks = function(y) c(min(y), max(y))) + 
  labs(x = "Easting (m)", y = "Northing (m)") +
  coord_equal()

p4 <- ggplot(chm_dsmtin_10_df, aes(x = x, y = y, fill = Z)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Triangulation over 10 m", fill = "Z") +
  scale_x_continuous(breaks = function(x) c(min(x), max(x))) +
  scale_y_continuous(breaks = function(y) c(min(y), max(y))) +   
  labs(x = "Easting (m)", y = "Northing (m)") +
  coord_equal()

p5 <- ggplot(chm_pitfree_df, aes(x = x, y = y, fill = Z)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Pit-free", fill = "Z") +
  scale_x_continuous(breaks = function(x) c(min(x), max(x))) +
  scale_y_continuous(breaks = function(y) c(min(y), max(y))) + 
  labs(x = "Easting (m)", y = "Northing (m)") +
  coord_equal() + 
  theme(legend.position = "right")

combined_plot <- wrap_plots(p1, p2, p5, ncol = 3) +
  plot_annotation(title = "Canopy Height Model (CHM)")  

# Print the combined plot
print(combined_plot)

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

# Locating trees
plot(chm_pitfree_subcirlce, col = magma(50))

## Cadaster
plot(sf::st_geometry(tree_cadaster),  add = TRUE, pch = 18, col="#F21B1B", cex = 3) #red

## LIDAR
ttops_LAS_ws_8 <- locate_trees(las_unfiltered, lmf(ws = 8, hmin=5))
plot(sf::st_geometry(ttops_LAS_ws_8), add = TRUE, pch = 1, col="#3647EB", cex = 2, lwd = 3) #blue
ttops_LAS_ws_deciduous <- locate_trees(las_unfiltered, lmf(ws_deciduous, hmin=5, shape = c("square"), ws_args = list("Z")))
plot(sf::st_geometry(ttops_LAS_ws_deciduous), add = TRUE, pch = 2, col="#00C87F", cex = 2, lwd = 3) #blue
ttops_LAS_ws_combined <- locate_trees(las_unfiltered, lmf(ws_combined, hmin=5,  ws_args = list("Z")))
plot(sf::st_geometry(ttops_LAS_ws_combined), add = TRUE, pch = 3, col="#ff0000", cex = 2, lwd = 3) #blue

## CHM
ttops_dsmtin_base <- locate_trees(chm_dsmtin_base, lmf(ws = 8, hmin=5))
plot(sf::st_geometry(ttops_dsmtin_base), add = TRUE, pch = 1, col="#00E3EC", cex = 2, lwd = 3) #lightblue
ttops_dsmtin_deciduous <- locate_trees(chm_dsmtin_base, lmf(ws_deciduous, hmin=5))
plot(sf::st_geometry(ttops_dsmtin_deciduous), add = TRUE, pch = 2, col="#00C87F", cex = 2, lwd = 3) #green
ttops_pitfree_subcirlce_deciduous <- locate_trees(chm_pitfree_subcirlce, lmf(ws = ws_deciduous, hmin = 5,  ws_args = list("Z")))
plot(sf::st_geometry(ttops_pitfree_subcirlce_deciduous), add = TRUE, pch = 3, col="#ff0000", cex = 2, lwd = 3) #green


## compared to the Cadaster, these 2 extra tree tops are not even trees
# CHM is faster because there are less data to process, 
# but it is also more complex because the output depends on how the CHM has been built.
## GGPLOT
'Color Theme Swatches in Hex
  .Primary-Colors-1-hex { color: #3647EB; }
      .Primary-Colors-2-hex { color: #00E3EC; }
          .Primary-Colors-3-hex { color: #F21B1B; } red
              .Primary-Colors-4-hex { color: #00C87F; }
                  .Primary-Colors-5-hex { color: #FFF200; } yellow'

chm_df <- as.data.frame(chm_pitfree_subcirlce, xy = TRUE)
names(chm_df)[3] <- "value"

# ✅ Get first and last coordinate values for axis breaks
x_breaks <- range(chm_df$x, na.rm = TRUE)
y_breaks <- range(chm_df$y, na.rm = TRUE)

# ✅ Convert spatial data to dataframe format (ensuring EPSG:25833)
convert_sf_to_df <- function(sf_obj) {
  sf_obj <- st_transform(sf_obj, crs = 25833)  # Ensure correct CRS
  df <- as.data.frame(st_coordinates(sf_obj))  # Extract x, y
  df <- df[complete.cases(df), ]  # Remove NAs
  return(df)
}

tree_cadaster_df <- convert_sf_to_df(tree_cadaster)
ttops_LAS_df <- convert_sf_to_df(ttops_LAS_ws_deciduous)
ttops_dsmtin_df <- convert_sf_to_df(ttops_dsmtin_deciduous)
ttops_pitfree_subcirlce_df <- convert_sf_to_df(ttops_pitfree_subcirlce_deciduous)

# ✅ Base ggplot layer (CHM + Cadaster Trees) in UTM EPSG:25833
base_map <- ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(option = "mako", name = "CHM") +
  geom_point(data = tree_cadaster_df, aes(x = X, y = Y), color = "#fc9a19", size = 4, shape = 18) +  # Convert to points
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_equal() +
  labs(x = "Easting (m)", y = "Northing (m)")

# 🟦 Compare cadaster vs LIDAR detected trees
plot_lidar <- base_map +
  geom_point(data = ttops_LAS_df, aes(x = X, y = Y), color = "#FFF200", size = 1, shape = 3,  stroke = 2) +
  labs(title = "Cadaster vs LiDAR Trees")

# 🔵 Compare cadaster vs CHM detected trees (DSM-TIN)
plot_chm_3 <- base_map +
  geom_point(data = ttops_dsmtin_df, aes(x = X, y = Y), color = "#FFF200", size = 1, shape = 3,  stroke = 2) +
  labs(title = "Cadaster vs CHM Trees (DSM-TIN 3)")

# 🟢 Compare cadaster vs CHM detected trees (pitfree)
plot_chm_5 <- base_map +
  geom_point(data = ttops_pitfree_subcirlce_df, aes(x = X, y = Y), color = "#FFF200", size = 1, shape = 3,  stroke = 2) +
  labs(title = "Cadaster vs CHM Trees (DSM-TIN 5)") +
  theme(legend.position = "right")

# 🎯 Arrange all plots with patchwork (3 columns layout)
combined_plot <- (plot_lidar | plot_chm_3 | plot_chm_5)

# ✅ Print the fixed combined plot
print(combined_plot)

# ---------------------------------------------
# Postprocessing CHM: try to remove traffic signs and lamps
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
plot(sf::st_geometry(ttops_pitfree_subcirlce_cleaned), add = TRUE, pch = 3,lwd = 3, col="#FFF200") #yellow

## GGPLOT
# ✅ Convert raster to dataframe for ggplot
raster_raw_df <- as.data.frame(rast(test), xy = TRUE)
names(raster_raw_df)[3] <- "value"  # Standardize raster value column name

# ✅ Convert `ttops_dsmtin_5_cleaned` to dataframe
ttops_pitfree_subcirlce_cleaned_df <- as.data.frame(st_coordinates(ttops_pitfree_subcirlce_cleaned))  # Extract x/y

# ✅ Get first and last coordinate values for axis breaks
x_breaks <- range(raster_raw_df$x, na.rm = TRUE)
y_breaks <- range(raster_raw_df$y, na.rm = TRUE)

# ✅ Create ggplot with Raw Raster + Points
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
             color = "#FFF200", size = 2, shape = 3, stroke = 2) +  # ✅ Yellow cross points
  scale_y_continuous(breaks = y_breaks) +
  labs(title = "Tree tops (postprocessed)", x = "Easting (m)", y = "Northing (m)") +
  coord_equal() 

# ✅ Print the plot
print(plot_raster_raw)
# ----------------------------------------------
# Get the stats of treetop detection
tree_cadaster <- tree_cadaster %>%
  mutate(buffer_size = ws_deciduous(baumhoehe))  # Apply function

# Create buffer zones
tree_cadaster_buffer <- st_buffer(tree_cadaster, dist = tree_cadaster$buffer_size)

plot(chm_pitfree_subcirlce_cleaned, col = mako(50))
plot(st_geometry(tree_cadaster_buffer), col = "transparent", border = "red", add = T, main = "Tree Buffers Based on Height")
plot(st_geometry(tree_cadaster), col = "red", pch = 16, add = TRUE)  # Add tree points

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


# ✅ Create the base map
base_map <- ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = tree_cadaster_buffer_poly, 
               aes(x = x, y = y, group = polygon_id), 
               color = "red", fill = NA, linewidth = 0.7) +
  geom_point(data = tree_cadaster_df, aes(x = X, y = Y), color = "red", size = 2) +
  scale_fill_viridis_c(option = "mako", name = "CHM") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  labs(x = "Easting (m)", y = "Northing (m)", title = "Tree Buffers Based on Height")

# ✅ Display the plot
print(base_map)
## -------------------------------

# Function to calculate true positives and false positives
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

# Calculate metrics for each method
results <- bind_rows(
  calculate_metrics(ttops_LAS_ws_8, "LIDAR ws=8"),
  calculate_metrics(ttops_LAS_ws_deciduous, "LIDAR ws=d"),
  calculate_metrics(ttops_dsmtin_deciduous, "CHM dsmtin ws=d"),
  calculate_metrics(ttops_pitfree_subcirlce_deciduous, "CHM pitfree ws=d"),
  calculate_metrics(ttops_pitfree_subcirlce_cleaned, "CHM pitfree ws=d [cleaned]")
)

# Sort results by True Positives in descending order
results_sorted <- results %>%
  arrange(desc(True_Positives))

# Convert Method to a factor with levels in sorted order
results_sorted$Method <- factor(results_sorted$Method, levels = results_sorted$Method)

# Reshape data for ggplot
results_long <- results_sorted %>%
  pivot_longer(cols = c(True_Positives, False_Positives), names_to = "Category", values_to = "Count")

# Plot Intersection Count with False Positives
ggplot(results_long, aes(x = Method, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Intersection Count with Cadaster",
       x = "Detection Method",
       y = "Count") +
  scale_fill_manual(values = c("True_Positives" = theme.main, "False_Positives" = theme.secondary)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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
ggplot(total_counts, aes(x = Method, y = Percentage_Deviation, fill = Percentage_Deviation > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Tree Count: Percentage Deviation from Cadaster",
       x = "Detection Method",
       y = "Percentage Deviation (%)",
       fill = "Deviation") +
  scale_fill_manual(values = c("TRUE" = theme.main, "FALSE" = theme.secondary)) +  # Green for overestimation, Red for underestimation
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
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
plot(sf::st_geometry(ttops), add = TRUE, pch = 3, col="#ff0000", cex = 2, lwd = 3) #green

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
ccm = ~custom_crown_metrics(z = Z, i = Intensity)

## Dalponte Method (raster-based)
algo_dalponte <- dalponte2016(chm, ttops,   th_tree = 2,    # Minimum tree height
                              th_seed = 0.45, # Seed threshold for initial growth
                              th_cr = 0.55,   # Crown merging threshold
                              max_cr = 20     # Max crown diameter (20 pixels = 10m for 0.5m CHM)
)
las_dalponte <- segment_trees(las_unfiltered, algo_dalponte)

#tree10 <- filter_poi(las_tree, treeID == 11)
#plot(tree10, size = 8, bg = "white")
crowns_dalponte <- crown_metrics(las_dalponte, func = ccm, geom = "concave")
## Post processing: Filter out small convex hulls, consider crowns over 1m2
plot(crowns_dalponte["z_max"], main = "x_max")
#crowns <- calculate_ratio_df(crowns)
#crowns <- crowns[crowns$width_to_height_ratio > 0.5 & crowns$width_to_height_ratio < 1.5 ,]
plot(chm, col = mako(50))
#plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
plot(crowns_dalponte["z_max"], col = "transparent", border = "#4baabf", lwd = 3, add = TRUE)
plot(sf::st_geometry(ttops), add = TRUE, pch = 3, col="#ff0000", cex = 2, lwd = 3)

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
plot(crowns_silva["i_mean"], col="transparent", border = "#9146fa", lwd = 3, add = TRUE)
plot(sf::st_geometry(ttops), add = TRUE, pch = 3, col="#ff0000",  cex = 2,lwd = 2)

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
  scale_fill_viridis_c(option = "mako", name = "CHM") +
  geom_point(data = tree_cadaster_df, aes(x = X, y = Y), color = "#fc9a19", size = 4, shape = 18) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_equal() +
  labs(x = "Easting (m)", y = "Northing (m)")

# Dalponte method visualization
plot_dalponte <- base_map +
  geom_polygon(data = crowns_dalponte_df, aes(x = x, y = y, group = polygon_id),
               color = "#4baabf", fill = NA, linewidth = 0.8) +
  geom_point(data = ttops_LAS_df, aes(x = X, y = Y), color = "red", shape = 3, size = 3) +
  labs(title = "Dalponte 2016 Segmentation")

# Silva method visualization
plot_silva <- base_map +
  geom_polygon(data = crowns_silva_df, aes(x = x, y = y, group = polygon_id),
               color = "#9146fa", fill = NA, linewidth = 0.8) +
  geom_point(data = ttops_LAS_df, aes(x = X, y = Y), color = "red", shape = 3, size = 3) +
  labs(title = "Silva 2016 Segmentation")

# Arrange both plots side by side
combined_plot <- (plot_dalponte | plot_silva)

# Print the fixed combined plot
print(combined_plot)

