library(lidR)
library(rgl)

las_nobuild_path <- "C:/tree-canopy/data/400_5816/LAS_no_buildings"
las_files <- list.files(path = las_nobuild_path, pattern = "\\.las$", full.names = TRUE, recursive = FALSE)

las <- readLAS(las_files[1])

plot(las, color = "Z", size = 1, bg = "white")
