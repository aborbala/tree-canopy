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
library(osmdata)

# Function to remove buildings from Lidar data
remove_structures_from_las  <- function(merged_structures_path, las_path) {
  structures <- st_read(merged_structures_path, quiet = TRUE)
  crs <- st_crs(structures)
  las <- readLAS(las_path, filter = "-drop_z_below 0")
  las_check(las)
  st_crs(las) <- 25833
  
  # Add 50 cm buffer around the buildings
  structures <- st_buffer(structures, 1)
  
  # Normalization (subtract the DTM)
  gnd <- filter_ground(las)
  dtm <- rasterize_terrain(las, 1, knnidw())
  nlas <- normalize_height(las, knnidw())
  
  # Get the inverse of the buildings and crop it
  las_bbox <- st_as_sfc(st_bbox(las))
  structures_aoi <- st_crop(structures, las_bbox)
  erase_feature <- st_difference(las_bbox, st_union(structures_aoi))

  las_clip <- clip_roi(nlas, erase_feature)
  
  return(las_clip)
}

get_bridges_osm <- function(buildings_path) {
  buildings <- st_read(buildings_path, quiet = TRUE)
  local_crs <- st_crs(buildings)
  
  # Convert bounding box to WGS 84 (EPSG:4326) for OSM query
  bbox_local <- st_bbox(buildings)
  bbox_wgs84 <- st_transform(st_as_sfc(bbox_local, crs = st_crs(buildings)), crs = 4326)
  bbox_wgs84 <- st_bbox(bbox_wgs84)  # Convert to bbox format
  
  # Query OSM for bridges within the bounding box
  query <- osmdata::opq(bbox = c(bbox_wgs84["xmin"], bbox_wgs84["ymin"], bbox_wgs84["xmax"], bbox_wgs84["ymax"])) %>%
    osmdata::add_osm_feature(key = "bridge", value = c("viaduct", "yes", "aqueduct", "beam", "suspension", "cantilever"))
  
  osm_bridges <- osmdata_sf(query)
  
  # Extract only the relevant geometries (lines/polygons)
  if (!is.null(osm_bridges$osm_lines)) {
    bridges_sf <- osm_bridges$osm_lines
  } else if (!is.null(osm_bridges$osm_polygons)) {
    bridges_sf <- osm_bridges$osm_polygons
  } else {
    warning("No bridges found in the bounding box.")
    return(NULL)
  }
  
  bridges_sf <- st_transform(bridges_sf, crs = local_crs)
  
  return(bridges_sf)
}

# Folder paths
las_folder_path <- "C:/tree-canopy/data/386_5818/LAS"
buildings_path <- "C:/tree-canopy/data/386_5818/buildings.gpkg"
bridges_output_path <- "C:/tree-canopy/data/386_5818/bridges.gpkg"
merged_structures_path <- "C:/tree-canopy/data/386_5818/structures.gpkg"
output_folder_path <- "C:/tree-canopy/data/386_5818/LAS_no_buildings"

if (!dir.exists(output_folder_path)) dir.create(output_folder_path)

bridges_sf <- get_bridges_osm(buildings_path)

if (!is.null(bridges_sf)) {
  # Overwrite existing bridges file
  st_write(bridges_sf, bridges_output_path, quiet = TRUE, delete_layer = TRUE)
  cat("Bridges extracted and saved to:", bridges_output_path, "\n")
  
  # Load buildings data
  buildings_sf <- st_read(buildings_path, quiet = TRUE)
  
  # Keep only geometries (drop attributes)
  buildings_geom <- st_geometry(buildings_sf)
  bridges_geom <- st_geometry(bridges_sf)
  
  # Union buildings and bridges into one geometry set
  merged_structures_geom <- st_union(buildings_geom, bridges_geom)
  
  # Convert back to sf object (without attributes)
  merged_structures <- st_sf(geometry = merged_structures_geom, crs = st_crs(buildings_sf))
  
  # Overwrite existing merged structures file
  st_write(merged_structures, merged_structures_path, quiet = TRUE, delete_layer = TRUE)
  cat("Merged structures saved to:", merged_structures_path, "\n")
} else {
  cat("No bridges found. Using only buildings dataset for removal.\n")
  merged_structures_path <- buildings_path  # Use buildings only if no bridges found
}


# Loop through all LAS files in the folder and remove buildings
las_files <- list.files(las_folder_path, pattern = "\\.las$", full.names = TRUE)
for (las_path in las_files) {
  cat("Processing:", las_path, "\n")
  las_nobuild <- remove_structures_from_las(merged_structures_path, las_path)
  
  # Save the resulting LAS file without buildings
  output_file_name <- paste0(file_path_sans_ext(basename(las_path)), "_nobuild.las")
  output_file_path <- file.path(output_folder_path, output_file_name)
  writeLAS(las_nobuild, output_file_path)
}

