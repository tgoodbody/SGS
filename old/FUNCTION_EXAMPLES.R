 sample_ahels(mraster = mr[[1:3]], 
              existing = e, 
              plot = TRUE)

 sample_ahels(mraster = mr[[1:3]], 
              existing = e, 
              nQuant = 20, 
              nSamp = 300,
              filename = tempfile(fileext = ".shp"))
 
 
 
calculate_coobs(mraster = mr,
                existing = e,
                cores = 4,
                details = TRUE,
                filename = tempfile(fileext = ".shp"))


 poplhs <- calculate_lhsPop(mraster = mr)

 calculate_lhsOpt(popLHS = poplhs)

 calculate_lhsOpt(popLHS = poplhs, 
                  PCA = FALSE, 
                  iter = 200)
 
  calculate_lhsPop(mraster = mr, 
                   nQuant = 10, 
                   PCA = FALSE)
  
  

   calculate_reqSamples(sraster = sr, 
                        nSamp = 200)
  
   e.sr <- extract_strata(sraster = sr, 
                          existing = e)
  
   calculate_reqSamples(sraster = sr, 
                        nSamp = 200, 
                        existing = e.sr)
   
   
    sample_clhs(mraster = mr, 
                nSamp = 200, 
                plot = TRUE, 
                iter = 100)
   
    sample_clhs(mraster = mr, 
                nSamp = 400, 
                existing = e, 
                iter = 250,
                details = TRUE)
   
    sample_clhs(mraster = mr, 
                nSamp = 200, 
                iter = 200, 
                existing = e,
                access = ac, 
                buff_inner = 100,
                buff_outer = 300, 
                plot = TRUE)
   
    #--- cost constrained examples ---#
    #--- calculate distance to access layer for each pixel in mr ---#
    mr.c <- calculate_distance(raster = mr, 
                               access = ac)
   
    sample_clhs(mraster = mr.c, 
                nSamp = 250, 
                iter = 200,
                cost = "dist2access",
                plot = TRUE)
   
    sample_clhs(mraster = mr.c, 
                nSamp = 250, 
                existing = e, 
                iter = 200,
                cost = "dist2access", 
                plot = TRUE)
    
    
     #--- perform grid sampling ---#
     sample_grid(raster = sr, 
                 gridsize = 200, 
                 plot = TRUE)
    
     sample_grid(raster = sr, 
                 gridsize = 100, 
                 access = ac,
                 buff_inner = 50,
                 buff_outer = 200)
    
     sample_grid(raster = sr,
                 gridsize = 200,
                 access = ac,
                 buff_inner = 100,
                 buff_outer = 400,
                 filename = tempfile(fileext = ".shp"),
                 plot = TRUE)
     
     
     
     
     
      #--- perform simple random sampling ---#
      sample_srs(raster = sr, 
                 nSamp = 200, 
                 plot = TRUE)
     
      sample_srs(raster = sr, 
                 nSamp = 200, 
                 access = ac,
                 mindist = 200,
                 buff_inner = 50,
                 buff_outer = 200)
     
      sample_srs(raster = sr,
                 nSamp = 200,
                 access = ac,
                 buff_inner = 50,
                 buff_outer = 200,
                 filename = tempfile(fileext = ".shp"))
      
      
      
       #--- perform stratified sampling random sampling ---#
       sample_strat(sraster = sr, 
                    nSamp = 200, 
                    plot = TRUE)
      
       sample_strat(sraster = sr, 
                    nSamp = 200, 
                    access = ac,
                    existing = e, 
                    mindist = 200, 
                    buff_inner = 50, 
                    buff_outer = 200)
      
       sample_strat(sraster = sr, 
                    nSamp = 200, 
                    access = ac,
                    buff_inner = 50, 
                    buff_outer = 200, 
                    filename = tempfile(fileext = ".shp"))
       
       
        strat_breaks(mraster = mr, 
                     metric = "zmax", 
                     breaks = br.max,
                     plot = TRUE, 
                     details = TRUE)
       
        strat_breaks(mraster = mr, 
                     metric = 1, 
                     metric2 = "zsd",
                     breaks = br.max, 
                     breaks2 = br.sd, 
                     plot = TRUE)
        
        
        
        
        
         #--- perform stratification using k-means ---#
         kmeans <- strat_kmeans(mraster = mr, 
                                nStrata = 5)
        
         kmeans <- strat_kmeans(mraster = mr, 
                                nStrata = 5, 
                                iter = 1000,
                                algorithm = "MacQueen",
                                plot = TRUE, 
                                details = TRUE)
        
         kmeans <- strat_kmeans(mraster = mr, 
                                nStrata = 5, 
                                iter = 1000,
                                plot = TRUE, 
                                filename = tempfile(fileext = ".tif"), 
                                overwrite = TRUE)
         
         
         
          #--- perform optimum sample boundary stratification ---#
          strat_osb(mraster = mr, 
                    metric = "zsd", 
                    nSamp = 200, 
                    nStrata = 4, 
                    plot = TRUE)
         
          strat_osb(mraster = mr, 
                    metric = 4, 
                    nSamp = 20, 
                    nStrata = 3, 
                    plot = TRUE, 
                    details = TRUE)
         
          strat_osb(mraster = mr, 
                    metric = "zmax", 
                    nSamp = 100, 
                    nStrata = 5, 
                    subset = 0.75, 
                    filename = tempfile(fileext = ".tif"))
          
          
           strat_pcomp(mraster = mr, 
                       nStrata = 5, 
                       plot = TRUE)
          
           strat_pcomp(mraster = mr, 
                       nStrata = 4, 
                       nStrata2 = 4, 
                       plot = TRUE, 
                       details = TRUE)
          
           strat_pcomp(mraster = mr, 
                       nStrata = 3, 
                       nStrata2 = 3, 
                       filename = tempfile(fileext = ".tif"))
           
           
           
           
            strat_quantiles(mraster = mr, 
                            metric = 4, 
                            nQuant = 10,
                            plot = TRUE, 
                            details = TRUE)
           
            strat_quantiles(mraster = mr, 
                            metric = "zsd",
                            metric2 = "zq95", 
                            nQuant = 3, 
                            nQuant2 = 4)
           
            strat_quantiles(mraster = mr, 
                            metric = 1, 
                            metric2 = "zsd",
                            nQuant = 2, 
                            nQuant2 = 2, 
                            filename = tempfile(fileext = ".tif"))