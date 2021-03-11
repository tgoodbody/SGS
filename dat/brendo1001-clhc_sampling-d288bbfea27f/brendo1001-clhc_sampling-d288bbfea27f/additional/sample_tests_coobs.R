## Algorithm for helping determine the locations of addtional samples given the existence of prior sampling. 
### This method is based on a count of obervations approach.
### We want to determine at every grid cell or pixel, how many observations are similar in terms of the covariate space.

###Notes
# The example below only considers numerical data 
# The example below has been created using some supplied data. Will need to be adapted for other datasets




#Libraries
library(raster);library(sp); library(rgdal); library(snow); library(doParallel); library(parallel)


#Given pre-exisiting samples how do we incorporate those and then do a cLHC sample?
setwd("Z:/Dropbox/2018/cLHC_samplingPAPER/gitHub/muddles/additional")

#Point data
dat<- read.table("HunterValley_SiteObsAll.txt", header = T,sep = ",")


#covariates
load(file = "HV_coobs.rda")
str(tempD)

#Make a raster stack of covariates
s1<- stack()
for (i in 4:ncol(tempD)){
  r1<- rasterFromXYZ(tempD[,c(1,2,i)])
  names(r1)<- names(tempD)[i]
  s1<- stack(s1,r1)}


#covariance matrix
covMat<- as.matrix(cov(tempD[,4:10]))

#extract covariate data at points
coordinates(dat)<- ~ X +Y 
DSM_data<- extract(s1,dat, sp= 1, method = "simple")
dat<- as.data.frame(DSM_data)
dat<- dat[complete.cases(dat),]

# Begin a parallel cluster and register it with foreach:
cpus = 2 # The number of nodes/cores to use in the cluster
cl <- makeCluster(spec = cpus)
# Register the cluster with foreach:
registerDoParallel(cl)


oper1 <- foreach(i=1:nrow(tempD), .combine = "c") %dopar% {  #para mode
#for (i in 1:nrow(tempD)){  #sequential mode
  
  # Pixel distance
  pix<- tempD[i,4:10] #pixel values
  pixDist<- mahalanobis(x = as.matrix(tempD[,4:10]), center = as.matrix(pix), cov =covMat) # calculate distance of selected pixel to all other pixels
  minPix<-min(pixDist) # minimum distance (will be 0 always)
  maxPix<- quantile(pixDist, probs = 0.975) #maximum distance  ##the probs variable could change   ####### Hack to avoid outliers
  
  #data distance
  datDist<- mahalanobis(x = as.matrix(dat[,6:12]), center = as.matrix(pix), cov =covMat) #calculate distance of observations to all other pixels
  datNdist<- (datDist-minPix)/(maxPix-minPix) # standardarise 
  
  datNdist[datNdist > 1] <- 1 #if the datNdist is greater than 1 that means it datDist is greater than maxDist ##HACK
  datNdist <- 1- datNdist  # Higher values mean more similar
  
  #count how many obs are above a given threshold
  sum(datNdist >= 0.975)}

#stop cluster 
stopCluster(cl)


#prepare grid outputs and export
tempD$sampleNOS<- oper1
r1<- rasterFromXYZ(tempD[,c("x", "y", "sampleNOS")])
plot(r1)
writeRaster(r1, filename="sampleNos.tif", format="GTiff", overwrite=TRUE)
  
  
