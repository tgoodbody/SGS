# This script provides an example of a method for relocating a sample point selected by cLHC sample to an alternative location
# that does not degrade the integrity of the sampling objective.
# Relocating a sample involves a few steps:
# 1. Set a buffer around the site that can not be visited from which potential alternative sampling locations can be assessed
# 2. Estimate the distance between selcted sites and all points within the buffer area. For numerical data this could be something like mahalanobis distance
# For categorical data we just remove the points in the buffer area that do not share the same categorical levels as the selcted site
# 3. Fit a membership function (this case a sigmoid curve) to the estimated distances
# 4. Select alternaitve sites where the membership is greater than or equal to a given threshold

### Note is the example below we use some example data. The script has not been generalised to the point where one could run the code without issue 
### given a different data set. Things like column indexes will need to be modified to run with another data set. The script below
### will work with the provided data and should be easy to adapt giving some other data.



library(clhs)
library(raster)
library(sp)
library(rasterVis)
library(entropy)

setwd("Z:/Dropbox/2017/pedometricsCOnference/clhs/cLHC_work/relocate") # set working directory accordingly


#load raster data
load(file="hvRasters.rda")
s1

################################################################################################
###Do a conditioned latin hypercube sample
#1. First convert rasters to data frame
tempD<- data.frame(cellNos = seq(1:ncell(s1)))
vals<- as.data.frame(getValues(s1))
tempD<- cbind(tempD, vals)
tempD<- tempD[complete.cases(tempD), ]
cellNos<- c(tempD$cellNos)
gXY<- data.frame(xyFromCell(s1, cellNos, spatial = F))
tempD<- cbind(gXY,tempD)
str(tempD)
#factor variables
tempD$dlwc_sl1<- as.factor(tempD$dlwc_sl1)
tempD$clippedGeol1<- as.factor(tempD$clippedGeol1)
str(tempD)

#2.Do cLHS
#res <- clhs(tempD[,4:10], size = 100, iter = 100000, progress = T, simple = TRUE)
#str(res)
################################################################################################
##Alternatively just load in a possible sampling configuration to save time
load(file="sampleTest.rda")
################################################################################################




################################################################################################
#Extract sample from dataframe
sampDat<- tempD[res,]
#plot
rx<- rasterFromXYZ(tempD[,c(1,2,9)])
plot(rx)
points(sampDat[,1:2], pch=20)
#remove sample points from selection frame
tempC<- tempD[-res,]
#########################################

## Coordinate reference systems
#coordinates(sampDat)<- ~ x + y 
#proj4string(sampDat) <- CRS("+init=epsg:32756")
#dat@proj4string

## Write point data to shapefile
#writeOGR(sampDat, ".", "HV_clhs_orig", "ESRI Shapefile")





#######
# Sigmoid function approach 
buff<- 500 #buffer radius set as metres
memo<- 0.99 # membership value (selection criteria)


#variance-covariance matrix (for mahalinobis distance calculation)
covM<-cov(tempD[,c(4,5,6,8,9)])
covM


##output for new re-located samples
outs1<- sampDat[1,]
outs1$dist<- NA
outs1$distances<- NA
outs1$zz<- NA
outs1<- outs1[0,]
outs1


#plot of area
plot(rx)


##########
## Beginning of selection of alterntive site
## Note in the example below, we loop through each site in turn and select an alternative site
## In a normal situation this would be really happen
## Instead selecting an alternative site would be on a site by site basis




# select a point at random 
for (i in 1:nrow(sampDat)) {
  points(sampDat[i,1:2], pch=20, cex=2) #row 1
  samp1<- sampDat[i,] #make object of selected site
  
  # Select all grid points within Xm of site
  dist<-spDistsN1(as.matrix(tempC[,1:2]), as.matrix(samp1[,1:2]), longlat = FALSE)  # the distance from function
  bind<-as.data.frame(cbind(tempC,dist))
  sorted <- subset(bind, dist < buff) # grid points within buffer zone

  ########
  #subset out those possibilties that do not have the same factor config as sample
  subs2<- which(sorted$dlwc_sl1==as.numeric(as.character(samp1$dlwc_sl1)) & sorted$clippedGeol1==as.numeric(as.character(samp1$clippedGeol1)))

  #extract those rows that fulfil categorical criteria
  sorted<- sorted[subs2,]
  
  # raster of the subset
  subR1<- rasterFromXYZ(sorted[,c(1:2,9)])
  plot(subR1)
  points(sampDat[i,1:2], pch=20, cex=2)
  
  ###### Estimate memberships
  #Mahalinobis distances of each grid cell to the selected location
  distances<- matrix(NA,ncol=1,nrow=nrow(sorted))  # matrix of where the distances go
  distances[,1]<- mahalanobis(as.matrix(sorted[,c(4,5,6,8,9)]), as.matrix(sampDat[i,c(4,5,6,8,9)]), cov= covM)
  
  #Append distances to possible sites
  sorted2<- cbind(sorted,distances)
  
  # raster of the subset
  subR2<- rasterFromXYZ(sorted2[,c(1:2,12)])
  plot(subR2) # plot the mahalanobis distance
  points(sampDat[i,1:2], pch=20, cex=2)
  
  #subset out possible sample locs 
  sorted3<- sorted2[order(sorted2$distances, decreasing = F),]  #sort by membership

  #Sigmoid function
  xx<- 1-(1/(1 + exp(-1*(sorted3$distances - median(sorted3$distances))))) # deviation from median
  zz<- (xx-min(xx))/(max(xx)-min(xx)) #normalise
  sorted3<- cbind(sorted3,zz)
  hist(sorted3$distances)
  plot(sorted3$distances,xx)
  plot(sorted3$distances,sorted3$zz, xlab="distance", ylab="similarity")
  
  # raster of the subset
  subR2<- rasterFromXYZ(sorted3[,c(1:2,13)])
  plot(subR2) # plot the membership with the buffer zone
  points(sampDat[i,1:2], pch=20, cex=2)
  
  #selection
  #1. Want all sites that fulfil membership criteria or at least the best 10 if there are less than 10 sites that fulfil criteria
  if (sum(sorted3$zz>=memo) <= 10) {sorted4<- sorted3[1:10, ]} else {sorted4<- sorted3[sorted3$zz>=memo,]}
  sorted4<- sorted4[order(sorted4$distances, decreasing = F),]  #sort by distance
  tr <- sample(nrow(sorted4), 1) # select a single site
  outs1[i,]<- sorted4[tr,]
  print (i)}
  




# Kullback-Leibler (KL)divergence

#First get the quantile matrix and bin count for the whole area

# covariate data
tempC_1<- tempC[,c(4,5,6,8,9)]
str(tempC_1)

# Number of bins
nosP<- 10

#quantile matrix (of the covariate data)
q.mat<- matrix(NA, nrow=(nosP+1), ncol= 5)
j=1
for (i in 1:ncol(tempC_1)){ #not the index start here
  #get a quantile matrix together of the covariates
  ran1<- max(tempC_1[,i]) - min(tempC_1[,i])
  step1<- ran1/nosP 
  q.mat[,j]<- seq(min(tempC_1[,i]), to = max(tempC_1[,i]), by =step1)
  j<- j+1}
q.mat




#covariate data hypercube
#############################################
## This takes a while to do so only do it once if you can 

cov.mat<- matrix(1, nrow=nosP, ncol=5)


for (i in 1:nrow(tempC_1)){ # the number of pixels
  cntj<- 1 
  for (j in 1:ncol(tempC_1)){ #for each column
    dd<- tempC_1[i,j]  
    for (k in 1:nosP){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}

cov.mat




####Compare whole study area covariate space with the orginal sample
dat<- sampDat
str(dat)
dat1<- dat[,c(4:6, 8:9)]
str(dat1)

#sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
h.mat<- matrix(1, nrow=nosP, ncol=5)

for (i in 1:nrow(dat1)){ # the number of observations
  cntj<- 1 
  for (j in 1:ncol(dat1)){ #for each column
    dd<- dat1[i,j]  
    for (k in 1:nosP){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat[k, cntj]<- h.mat[k, cntj] + 1}
    }
    cntj<- cntj+1
  }
}

h.mat 

#Kullback-Leibler (KL) divergence
klo.1<- KL.empirical(c(cov.mat[,1]), c(h.mat[,1])) #twi
klo.1
klo.2<- KL.empirical(c(cov.mat[,2]), c(h.mat[,2])) # slope
klo.2
klo.3<- KL.empirical(c(cov.mat[,3]), c(h.mat[,3])) #mrvbf
klo.3
klo.4<- KL.empirical(c(cov.mat[,4]), c(h.mat[,4])) # directInSol
klo.4
klo.5<- KL.empirical(c(cov.mat[,5]), c(h.mat[,5])) #DEM
klo.5
klo<- mean(c(klo.1, klo.2,klo.3,klo.4,klo.5))
klo  # value of 0 means no divergence





####Compare whole study area covariate space with the adjusted sample
dat<- outs1
str(dat)
dat2<- dat[,c(4:6, 8:9)]
str(dat2)

#sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
j.mat<- matrix(1, nrow=nosP, ncol=5)

for (i in 1:nrow(dat2)){ # the number of observations
  cntj<- 1 
  for (j in 1:ncol(dat2)){ #for each column
    dd<- dat2[i,j]  
    for (k in 1:nosP){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){j.mat[k, cntj]<- j.mat[k, cntj] + 1}
    }
    cntj<- cntj+1
  }
}

j.mat 

#Kullback-Leibler (KL) divergence
kla.1<- KL.empirical(c(cov.mat[,1]), c(j.mat[,1])) #twi
kla.1
kla.2<- KL.empirical(c(cov.mat[,2]), c(j.mat[,2])) # slope
kla.2
kla.3<- KL.empirical(c(cov.mat[,3]), c(j.mat[,3])) #mrvbf
kla.3
kla.4<- KL.empirical(c(cov.mat[,4]), c(j.mat[,4])) # directInSol
kla.4
kla.5<- KL.empirical(c(cov.mat[,5]), c(j.mat[,5])) #DEM
kla.5
kla<- mean(c(kla.1, kla.2,kla.3,kla.4,kla.5))
kla # value of 0 means no divergence



####Compare original sample with adjusted sample
#Kullback-Leibler (KL) divergence
kloa.1<- KL.empirical(c(h.mat[,1]), c(j.mat[,1])) #twi
kloa.1
kloa.2<- KL.empirical(c(h.mat[,2]), c(j.mat[,2])) # slope
kloa.2
kloa.3<- KL.empirical(c(h.mat[,3]), c(j.mat[,3])) #mrvbf
kloa.3
kloa.4<- KL.empirical(c(h.mat[,4]), c(j.mat[,4])) # directInSol
kloa.4
kloa.5<- KL.empirical(c(h.mat[,5]), c(j.mat[,5])) #DEM
kloa.5
kloa<- mean(kloa.1, kloa.2,kloa.3,kloa.4,kloa.5)
kloa



#histograms
#twi
hist(tempC_1$twi, breaks = c(q.mat[,1]), freq = TRUE, main = "twi: whole area")
hist(dat1$twi, breaks = c(q.mat[,1]), freq = TRUE, main = "twi: original sample")
hist(dat2$twi, breaks = c(q.mat[,1]), freq = TRUE, main = "twi: relocated sample")

#slope
hist(tempC_1$slope, breaks = c(q.mat[,2]), freq = TRUE, main = "slope: whole area")
hist(dat1$slope, breaks = c(q.mat[,2]), freq = TRUE, main = "slope: original sample")
hist(dat2$slope, breaks = c(q.mat[,2]), freq = TRUE, main = "slope: relocated sample")

#mrvbf
hist(tempC_1$mrvbf, breaks = c(q.mat[,3]), freq = TRUE, main = "mrvbf: whole area")
hist(dat1$mrvbf, breaks = c(q.mat[,3]), freq = TRUE, main = "mrvbf: original sample")
hist(dat2$mrvbf, breaks = c(q.mat[,3]), freq = TRUE, main = "mrvbf: relocated sample")

#insolation
hist(tempC_1$directInsol, breaks = c(q.mat[,4]), freq = TRUE, main = "insol: whole area")
hist(dat1$directInsol, breaks = c(q.mat[,4]), freq = TRUE, main = "insol: original sample")
hist(dat2$directInsol, breaks = c(q.mat[,4]), freq = TRUE, main = "insol: relocated sample")

#DEM
hist(tempC_1$DEM_streamBurned, breaks = c(q.mat[,5]), freq = TRUE, main = "elevation: whole area")
hist(dat1$DEM_streamBurned, breaks = c(q.mat[,5]), freq = TRUE, main = "elevation: original sample")
hist(dat2$DEM_streamBurned, breaks = c(q.mat[,5]), freq = TRUE, main = "elevation: relocated sample")


#save.image("clhs_relocate.RData") #save R session




