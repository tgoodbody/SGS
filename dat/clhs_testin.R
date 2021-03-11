https://www.sciencedirect.com/science/article/pii/S001670612030135X#f0005
library(spsann)


d <- df %>% sample_n(10,000,replace=F)

schedule <- scheduleSPSANN(chains = 1,
                           chain.length = 4,
                           initial.acceptance = 0.9,
                           x.max = 100,
                           y.max = 100,
                           temperature.decrease = 0.95)

set.seed(2001)
res <- optimCLHS(points = 100, 
                 candi = d[1:2], 
                 covars = d[3], 
                 use.coords = TRUE,
                 clhs.version = "update", 
                 weights = list(O1 = 0.5, O3 = 0.5), 
                 schedule = schedule)

objSPSANN(res) - objCLHS(points = res, 
                         candi = d[1:2], 
                         covars = d[3], 
                         use.coords = TRUE,
                         clhs.version = "paper", 
                         weights = list(O1 = 0.5, O3 = 0.5))

plot(res)

wall_poly_roads <- readRDS("wall_poly_roads.rds")

plot(wall_poly_roads$avg)

coords <- st_as_sf(res$points, coords = c("x","y"))

plot(coords,add=T)


library(LICORS)

kk <- kmeanspp(d,nstart = 100,start = "random",k=5)
kk

plot(df, col= kk$cluster + 1)
