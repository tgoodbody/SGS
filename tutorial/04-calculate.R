#################################
########## CALCULATE ############
#################################

##### Calculate distance from road access #####

#--- valuable as a cost variable for some sampling schemes ---#
calculate_distance(raster = sraster, # input
                   access = access, # define access road network
                   plot = TRUE) # plot

calculate_distance(raster = mraster, # input
                   access = access, # define access road network
                   filename = tempfile(fileext = ".tif")) # write file to disc

##### Calculate principal components #####

calculate_pcomp(mraster = mraster, # input
                nComp = 5, # number of components to output
                plot = TRUE ) # plot

calculate_pcomp(mraster = mraster, # input
                nComp = 3, # number of components to output
                plot = TRUE, # plot
                details = TRUE) # details about the principal component analysis appended

##### Calculate number of required samples #####
# right now includes proportional allocation - optimal allocation to be added

#--- perform grid sampling ---#
calculate_reqSamples(sraster = sraster, # input
                     nSamp = 200) # number of desired samples

calculate_reqSamples(sraster = sraster, # input
                     nSamp = 200) # number of desired samples

#--- calculate existing samples to include ---#
e.sr <- extract_strata(sraster = sraster, # input
                       existing = existing) # existing samples

calculate_reqSamples(sraster = sraster, # input
                     nSamp = 200, # number of desired samples
                     existing = e.sr) # existing samples

calculate_reqSamples(sraster = sraster, # input
                     nSamp = 200, # number of desired samples
                     existing = e.sr, # existing samples
                     include = TRUE) # include existing samples within total




sample_ahels(mraster = mraster[[1:3]], # input mraster - first 3 layers only
             existing = existing, # existing samples
             plot = TRUE) # plot



sample_ahels(mraster = mraster[[1:3]], # input mraster - first 3 layers only
             existing = existing, # existing samples
             nQuant = 20, # define 20 quantiles
             nSamp = 300, # total samples desired
             filename = tempfile(fileext = ".shp")) # write samples to disc


### A few more algorithms in the calculate_* category but they take a bit longer to run so ill show them another time.

#################################
########## CALCULATE ############
#################################
