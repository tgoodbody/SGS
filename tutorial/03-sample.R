#################################
########### SAMPLING ############
#################################

##### simple random sampling #####

#--- perform simple random sampling ---#
sample_srs(raster = sraster, # input sraster
           nSamp = 200, # number of desired samples
           plot = TRUE) # plot


#--- introduce access, min distance, and buffers ---#
sample_srs(raster = mraster, # input mraster
           nSamp = 200, # number of desired samples
           access = access, # define access road network
           mindist = 200, # minimum distance samples must be apart from one another
           buff_inner = 50, # inner buffer - no samples within this distance from road
           buff_outer = 200, # outer buffer - no samples further than this distance from road
           plot = TRUE) # plot


sample_srs(raster = sraster, # input
           nSamp = 200, # number of desired samples
           access = access, # define access road network
           buff_inner = 100, # inner buffer - no samples within this distance from road
           buff_outer = 400, # outer buffer - no samples further than this distance from road
           plot = TRUE, # plot
           filename = tempfile(fileext = ".shp")) # write output samples to file

##### grid sampling #####

#--- perform grid sampling ---#
sample_systematic(raster = sraster, # input sraster
            cellsize = 200, # grid distance
            plot = TRUE) # plot

sample_systematic(raster = sraster, # input sraster
            cellsize = 100, # grid distance
            access = access, # define access road network
            buff_inner = 50, # inner buffer - no samples within this distance from road
            buff_outer = 200) # outer buffer - no samples further than this distance from road

sample_systematic(raster = mraster, # input mraster
            cellsize = 200, # grid distance
            access = access, # define access road network
            buff_inner = 100, # inner buffer - no samples within this distance from road
            buff_outer = 400, # outer buffer - no samples further than this distance from road
            filename = tempfile(fileext = ".shp"), # write output samples to file
            plot = TRUE) # plot


##### stratified sampling #####

#--- perform stratified sampling random sampling ---#
sample_strat(sraster = sraster, # input sraster
             nSamp = 200, # desired sample number
             plot = TRUE) # plot


#--- an 'existing' plot network can be included ---#

#--- extract strata values to existing samples ---#              
e.sr <- extract_strata(sraster = sraster, # input sraster
                       existing = existing) # existing samples to add strata value to

e.sr


sample_strat(sraster = sraster, # input sraster
             nSamp = 200, # desired sample number
             access = access, # define access road network
             existing = e.sr, # existing samples with strata values
             mindist = 200, # minimum distance samples must be apart from one another
             buff_inner = 50, # inner buffer - no samples within this distance from road
             buff_outer = 200, # outer buffer - no samples further than this distance from road
             plot = TRUE) # plot



sample_strat(sraster = sraster, # input
             nSamp = 200, # desired sample number
             access = access, # define access road network
             existing = e.sr, # existing samples with strata values
             include = TRUE, # include existing plots in nSamp total
             buff_inner = 50, # inner buffer - no samples within this distance from road
             buff_outer = 200, # outer buffer - no samples further than this distance from road
             filename = tempfile(fileext = ".shp"), # write output samples to file
             plot = TRUE) # plot

##### conditioned latin hypercube sampling #####

sample_clhs(mraster = mraster, # input
            nSamp = 200, # desired sample number
            plot = TRUE, # plot 
            iter = 100) # number of iterations


sample_clhs(mraster = mraster, # input
            nSamp = 400, # desired sample number
            existing = existing, # existing samples
            iter = 100, # number of iterations
            details = TRUE) # clhs details


sample_clhs(mraster = mraster, # input
            nSamp = 300, # desired sample number
            iter = 100, # number of iterations
            existing = existing, # existing samples
            access = access, # define access road network
            buff_inner = 100, # inner buffer - no samples within this distance from road
            buff_outer = 300, # outer buffer - no samples further than this distance from road
            plot = TRUE) # plot



#--- cost constrained examples ---#
#--- calculate distance to access layer for each pixel in mr ---#
mr.c <- calculate_distance(raster = mraster, # input
                           access = access,
                           plot = TRUE) # define access road network

sample_clhs(mraster = mr.c, # input
            nSamp = 250, # desired sample number
            iter = 100, # number of iterations
            cost = "dist2access", # cost parameter - name defined in calculate_distance()
            plot = TRUE) # plot


sample_clhs(mraster = mr.c, # input
            nSamp = 250, # desired sample number
            existing = existing, # existing samples
            iter = 100, # number of iterations
            cost = "dist2access", # cost parameter - name defined in calculate_distance()
            plot = TRUE) # plot


##### balanced sampling #####

sample_balanced(mraster = mraster, # input
                nSamp = 200, # desired sample number
                plot = TRUE) # plot

sample_balanced(mraster = mraster, # input
                nSamp = 100, # desired sample number
                algorithm = "lcube", # algorithm type
                access = access, # define access road network
                buff_inner = 50, # inner buffer - no samples within this distance from road
                buff_outer = 200) # outer buffer - no samples further than this distance from road

#################################
########### SAMPLING ############
#################################
