#################################
########### STRATIFY ############
#################################

##### strat_kmeans #####

#--- perform stratification using k-means ---#
strat_kmeans(mraster = mraster, # input
             nStrata = 5) # algorithm will produce 4 strata

#--- introduce plot and details ---#
strat_kmeans(mraster = mraster, # input
             nStrata = 10, # algorithm will produce 10 strata
             iter = 1000, # set minimum number of interations to determine kmeans centers
             algorithm = "MacQueen", # use MacQueen algorithm
             plot = TRUE, # plot output
             details = TRUE) # output details - kmeans stratification data and output sraster

#--- introduce filename and overwrite ---#
strat_kmeans(mraster = mraster, # input
             nStrata = 5, # algorithm will produce 4 strata
             center = FALSE, # do not center data
             scale = FALSE, # do not scale data
             plot = TRUE, # plot output
             filename = tempfile(fileext = ".tif"), # write output sraster to file
             overwrite = TRUE) # overwrite file on disc if it exists

##### strat_pcomp #####

#--- perform stratification using principal components ---#
strat_pcomp(mraster = mraster, # input
            nStrata = 5, # 5 strata with primary PC only
            plot = TRUE) # plot

#--- introduce bivariate stratification ---#
strat_pcomp(mraster = mraster, # input
            nStrata = 4, # 4 strata with primary
            nStrata2 = 4, # 4 strata with secondary PC - will produce 16 output strata
            plot = TRUE, # plot
            details = TRUE) # produce output details

strat_pcomp(mraster = mraster, # input
            nStrata = 3, # 3 strata with primary PC
            nStrata2 = 3, # 4 strata with secondary PC - will produce 9 output strata
            filename = tempfile(fileext = ".tif")) # write output sraster to file

##### strat_pcomp #####

#--- define vector breaks ---#
br.max <- c(3, 5, 11, 18) #zmax breaks
br.sd <- c(1, 2, 5) #zsd breaks

#--- perform stratification using breaks ---#
strat_breaks(mraster = mraster, # input
             metric = "zmax", # covariate of interest - numeric index or character string
             breaks = br.max, # breaks for primary covariate
             plot = TRUE, # plot
             details = TRUE) # output details


strat_breaks(mraster = mraster, # input
             metric = 1, # primary covariate of interest - numeric index or character string
             metric2 = "zsd", # secondary covariate of interest - numeric index or character string
             breaks = br.max, # breaks for primary covariate
             breaks2 = br.sd, # breaks for secondary covariate
             plot = TRUE) # plot


##### strat_quantiles #####

#--- perform stratification using quantiles ---#
strat_quantiles(mraster = mraster, # input
			          metric = 4, # primary covariate of interest - numeric index or character string
			          nQuant = 10, # number of quantiles for primary covariate
			          plot = TRUE) # plot

strat_quantiles(mraster = mraster, # input
			          metric = "zsd", # primary covariate of interest - numeric index or character string
			          metric2 = "zq95", # secondary covariate of interest - numeric index or character string
			          nQuant = 3, # number of quantiles for primary covariate
			          nQuant2 = 4, # number of quantiles for secondary covariate
			          plot = TRUE) # plot

strat_quantiles(mraster = mraster, # input
			          metric = 1, # primary covariate of interest - numeric index or character string
			          metric2 = "zsd", # secondary covariate of interest - numeric index or character string
			          nQuant = 2, # number of quantiles for primary covariate
			          nQuant2 = 2, # number of quantiles for secondary covariate
			          filename = tempfile(fileext = ".tif")) # write output sraster to file


##### strat_osb #####

#Takes a bit longer to run... so I don't show in this demo :)

