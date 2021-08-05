#################################
#### PREPARE EXAMPLE DATASET ####
#################################

#--- Load neccecary libraries ---#

library(sgsR)
library(terra)
library(sf)

#--- set seed ---#
set.seed(2021)

#--- Load mraster and access files ---#
r <- system.file("extdata", "wall_metrics_small.tif", package = "sgsR")

#--- load the mraster using the terra package ---#
mraster <- terra::rast(r)
terra::plot(mraster)

#--- Access shapefiles ---#
a <- system.file("extdata", "roads.shp", package = "sgsR")

#--- load the access vector using the sf package ---#
access <- sf::st_read(a)

terra::plot(mraster[[1]])
terra::plot(access, add = TRUE, col = "black")


#### prepare a stratified raster example ####

#--- apply kmeans algorithm to metrics raster ---#
sraster <- strat_kmeans(mraster = mraster, # use mraster as input for sampling
                        nStrata = 4, # algorithm will produce 4 strata
                        plot = TRUE) # algorithm will plot output

#### prepare a existing sample dataset ####

#--- apply kmeans algorithm to metrics raster ---#
existing <- sample_srs(raster = mraster, # use mraster as input for sampling
                       nSamp = 200, # request 200 samples be taken
                       mindist = 100, # define that samples must be 100 m apart
                       plot = TRUE) # algorithm will plot output


#################################
#### PREPARE EXAMPLE DATASET ####
#################################

