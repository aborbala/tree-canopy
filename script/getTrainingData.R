library("rstudioapi")
library("Rcpp")
library("sf")
library("lidR")
library("stars")
library("terra")
library("raster")

# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
wd <- getwd()

#########

### Remove buildings from Lidar data ###

# Load Buildings and Lidar data
buildings <- st_read("../data/test 382_5826/buildings.gpkg")
crs <- st_crs(buildings)
las <- readLAS("../data/test 382_5826/3dm_33_382_5826_1_be.las", filter = "-drop_z_below 0")
las_check(las)
st_crs(las) <-25833

## Normalization (substract the DTM)
gnd <- filter_ground(las)
#plot(gnd, size = 3, bg = "white", color = "Classification")
dtm <- rasterize_terrain(las, 1, knnidw())
#plot(dtm, col = gray(1:50/50))
#nlas <- normalize_height(las, tin(), dtm = dtm)
nlas <- normalize_height(las, knnidw())
hist(filter_ground(nlas)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")#nlas <- las - dtm
plot(nlas)

## Inverse of the buildings:
las_bbox <- st_as_sfc(st_bbox(las))
#plot(buildings[2], add=T)

#ras <- st_as_sf(read_stars("../data/test 382_5826/dop20rgbi_33_382_5826_2_be_2022.jp2"))

#plot(buildings[1], add=T)
buildings_aoi <- st_crop(buildings, las_bbox)
#plot(buildings_aoi)
erase_feature <- st_difference(las_bbox, st_union(buildings_aoi))
#plot(erase_feature, col='salmon')
las_clip <- clip_roi(nlas, erase_feature)
#writeLAS(las_clip,"../data/test 382_5826/las_clipped.las")

################################################################################


# Khosravipour et al. pitfree algorithm ## TODO: read about this
thr <- c(0,2,5,10,15)
edg <- c(0, 1.5)

## Lidar without buildings
las_nobuild <- readLAS("../data/test 382_5826/las_clipped.las")
las_clip <- filter_poi(las_nobuild, Z >= 10)
chm <- rasterize_canopy(las_clip, 1, pitfree(thr, edg))
plot(chm)
getCrowns <- function (las_clip, bbox){
  # Get aoi, clip the las based on the raster size(bbox)
  las_aoi <- clip_roi(las_clip, bbox)
  #st_crs(las_aoi) <- 25833#
  if (dim(las_aoi@data)[1] < 5){
    return()
  }

  chm <- rasterize_canopy(las_aoi, 1, pitfree(thr, edg))
  #plot(chm)
  #plot(buildings_aoi[1], add=T)
  
  #########################
  # Detection
  #ttops <- locate_trees(las_clip, lmf(ws = 5))
  
  #plot(chm, col = height.colors(50))
  #plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
  
  chm_p2r_05 <- rasterize_canopy(las_aoi, 0.5, p2r(subcircle = 0.2), pkg = "terra")
  kernel <- matrix(1,3,3)
  chm_p2r_05_smoothed <- terra::focal(chm_p2r_05, w = kernel, fun = median, na.rm = TRUE)
  ttops_chm_p2r_05 <- locate_trees(chm_p2r_05, lmf(5))
  ttops_chm_p2r_05_smoothed <- locate_trees(chm_p2r_05_smoothed, lmf(5))
  
  par(mfrow=c(1,2))
  col <- height.colors(50)
  #plot(chm_p2r_05, main = "CHM P2R 0.5", col = col); plot(sf::st_geometry(ttops_chm_p2r_05), add = T, pch =3)
  #plot(chm_p2r_05_smoothed, main = "CHM P2R 0.5 smoothed", col = col); plot(sf::st_geometry(ttops_chm_p2r_05_smoothed), add = T, pch =3)
  
  
  # Segmentation
  algo <- dalponte2016(chm_p2r_05_smoothed, ttops_chm_p2r_05_smoothed)
  las_tree <- segment_trees(las_aoi, algo) # segment point cloud
  #plot(las_tree, bg = "white", size = 4, color = "treeID") # visualize trees
  ## Convex hulls
  crowns <- crown_metrics(las_tree, func = .stdtreemetrics, geom = "convex")
  #plot(crowns["convhull_area"], main = "Crown area (convex hull)")
  return(crowns)
}

# Get crowns for the whole area
#crowns <- getCrowns(nlas, erase_feature)
#st_write(crowns, "../data/test 382_5826/crowns.gpkg")

# Get crowns for small tif tiles
files <- list.files(path="G:/Meine Ablage/data/382_5826_1/sliced_imgs/", pattern="*.tif", full.names=TRUE, recursive=FALSE)
df <- data.frame((matrix(ncol = 4, nrow = 0)))
colnames(df) <- c("xmin", "xmax", "ymin", "ymax")
# df <- lapply(files, function(x) {
#   ras <- raster::raster(x) # load file
#   extent <- raster::extent(ras)
#   # lidR: xleft, ybottom, xright
#   ext <- c(extent@xmin, extent@xmax, extent@ymin, extent@ymax)
#   return(ext)
# })

lapply(files, function(x) {
  print(x)
  ras <- st_as_sf(read_stars(x))
  st_crs(ras) <- 25833
  ext <- st_bbox(ras)
  crowns <- getCrowns(las_clip, ext)
  if (is.null(crowns) == F){
    filename <- strsplit(sub(".*/", "", x), "\\.")[[1]][1]
    st_write(crowns, paste("../data/test 382_5826/382_5826_1_crowns/",filename,".geojson", sep=""), delete_dsn = T)
  }
  gc()
})


############### CUSTOMIZE TREE DETECTION
test <- "G:/Meine Ablage/data/382_5826_1/sliced_imgs/382_5826_1_6.tif"
ras_t <- st_as_sf(read_stars(test))
ext <- st_bbox(ras_t)
start <- Sys.time()
las_aoi <- clip_roi(las_clip, ext)
end <- Sys.time()
duration <- end - start
duration
chm <- rasterize_canopy(las_aoi, 1, pitfree(thr, edg))
plot(chm)

chm_p2r_05 <- rasterize_canopy(las_aoi, 0.5, p2r(subcircle = 0.2), pkg = "terra")
plot(chm_p2r_05)
kernel <- matrix(1,3,3)
chm_p2r_05_smoothed <- terra::focal(chm_p2r_05, w = kernel, fun = median, na.rm = TRUE)
ttops_chm_p2r_05 <- locate_trees(chm_p2r_05, lmf(5))
plot(ttops_chm_p2r_05)
ttops_chm_p2r_05_smoothed <- locate_trees(chm_p2r_05_smoothed, lmf(5))

par(mfrow=c(1,2))
col <- height.colors(50)
plot(chm_p2r_05, main = "CHM P2R 0.5", col = col); plot(sf::st_geometry(ttops_chm_p2r_05), add = T, pch =3)
plot(chm_p2r_05_smoothed, main = "CHM P2R 0.5 smoothed", col = col); plot(sf::st_geometry(ttops_chm_p2r_05_smoothed), add = T, pch =3)

# Segmentation
algo <- dalponte2016(chm_p2r_05_smoothed, ttops_chm_p2r_05_smoothed)
las_tree <- segment_trees(las_aoi, algo) # segment point cloud
plot(las_tree, bg = "white", size = 4, color = "treeID") # visualize trees
## Convex hulls
crowns <- crown_metrics(las_tree, func = .stdtreemetrics, geom = "convex")
plot(crowns)
#st_write(crowns, "../data/test 382_5826/crowns_test.geojson", delete_dsn = T)
#st_crs(las_aoi) <- 25833

# ras <- raster::raster(x) # load file
# extent <- raster::extent(ras)
# ext <- c(extent@xmin, extent@xmax, extent@ymin, extent@ymax)

