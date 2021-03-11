# Algorithm for determining how to allocate additional samples within a study area given some existing samples
# The method is based mainly on the Carre et al (2007) HELS and  algorithm
# Hypercube Evaluation of a Legacy Sample (HELS)
#
# basically it entails:
# 1. From the covariate data generate a matrix of quantiles (n x p) n = number of samples that are needed. p = number of available covariates
# 2. Using the quantile matrix create Hypercube matrices for both the existing point data and covariates.
# 3. Work out the densities of the elements in the hypercubes (count of obs/pixels within each quantile / total number of data or pixels)
# 4. Work out the ratio of sampling and grid densities... This will identify under and over sampling in the hypercube
# 5. To add additional samples:
#         1. rank the ratios from smallest to largest
#         2. workout the number of samples required to equalise quantile density of grids and sample data
#         3. Repeat step 5.2 until total number of additonal samples have been allocated.
## Created: 18/05/17

#load r session with all the necesary objects
# load("clhs_samp.RData")


# Libraries
library(raster);library(sp); library(rgdal); library(tidyverse)


#Given pre-exisiting samples how do we incorporate those and then do a cLHC sample?
setwd("G:/Documents/post_doc/SGS/dat/brendo1001-clhc_sampling-d288bbfea27f/brendo1001-clhc_sampling-d288bbfea27f/additional")


# Number of addtional samples to take
nosP<- 15

tempD <- rast_dt %>% dplyr::select(x,y,avg,std,cov,p99) %>% na.omit() %>% as.data.frame()



##INPUT DATA
#########################################################################################################################
#covariates
str(tempD)
names(tempD)



#quantile matrix (of the covariate data)
q.mat<- matrix(NA, nrow=(nosP+1), ncol= 4)
j=1
for (i in 3:ncol(tempD)){ #not the index start here
  #get a quantile matrix together of the covariates
  ran1<- max(tempD[,i]) - min(tempD[,i])
  step1<- ran1/nosP 
  q.mat[,j]<- seq(min(tempD[,i]), to = max(tempD[,i]), by =step1)
  j<- j+1
  }
q.mat




#covariate data hypercube
#############################################
## This takes a while to do so only do it once if you can 

cov.mat<- matrix(0, nrow=nosP, ncol=4)


for (i in 1:nrow(tempD)){ # the number of pixels
  cntj<- 1 
  for (j in 3:ncol(tempD)){ #for each column
    dd<- tempD[i,j]  
    for (k in 1:nosP){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}

cov.mat[which(cov.mat==0)]<- 0.000001 # small number so we dont have to deal with zeros
cov.mat # hypercube matrix
covDens<- cov.mat/nrow(tempD)
covDens # density matrix


#####################################################################################################################


dat <- shapefile("G:\\Documents\\post_doc\\RMFinventory\\inst\\extdata\\existing_plots.shp")

dat <- as(dat,"data.frame")
names(dat) <- c("plotID","x","y")


#Point data
# dat<- read.table("HunterValley_SiteObsAll.txt", header = T,sep = ",")  # existing soil point data

wall_metrics <- brick(system.file("extdata","wall_metrics_small.tif", package = "RMFinventory"))
names(wall_metrics) <- c("avg", "cov", "std","p10", "p20","p50","p70","p95", "p99","d0","d2","d4","dns")

#extract covariate data at points
coordinates(dat)<- ~ x +y 
DSM_data<- raster::extract(wall_metrics,dat, sp= 1, method = "simple") #extract
dat<- as.data.frame(DSM_data) %>% dplyr::select(x,y,avg,std,cov,p99)
dat<- dat[complete.cases(dat),]
# dat <- dat %>% dplyr::select(x,y,avg,std,cov,p99)


#Point dat (little data manipulations)
dat
names(dat)[1:2]<- c("x", "y")
# dat$cellNos<- 0
dat
nrow(dat)
ncol(dat)


#sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
h.mat<- matrix(0, nrow=nosP, ncol=4)

for (i in 1:nrow(dat)){ # the number of observations
  cntj<- 1 
  for (j in 3:ncol(dat)){ #for each column
    dd<- dat[i,j]  
    for (k in 1:nosP){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat[k, cntj]<- h.mat[k, cntj] + 1}
    }
    cntj<- cntj+1
  }
}



#density
h.mat[which(h.mat==0)]<- 0.000001  # small number so we dont have to deal with zeros
h.mat 
datDens<-h.mat/nrow(dat) # data density



#### selecting new samples

rat<- datDens/covDens # ratio of data density and covariate density
or<- order(rat) # rank the index where the biggest discrepancy is
or

## indexes of quantiles that are not adequately sampled ie where rat is less than 1
l1<- which(rat < 1, arr.ind = T)
l1<- cbind(l1,which(rat < 1) )
l1
length(rat)
nrow(l1)


# What is the level of the greatest discrepancy? (This is important for comparing to later on when we have the additional sample)
indy<- which(l1[,3]==or[1])
rp<- l1[indy, 1]
rc<- l1[indy, 2]
rat[rp, rc]
h.mat[rp, rc]
cov.mat[rp, rc] # nuber of pixel (biggest discrepancy)
datDens[rp, rc] # data density
covDens[rp, rc] # covariate density




#start from the highest discrepancy then work our way down
upSamp<- nosP
rpos<- 1
base<- 1 # constant so that covariate columns are easily to select (realted to the column positions of the covariates in the tempD data frame)


while (upSamp != 0){  # while the number of samples to allocate is greater than 0
  indy<- which(l1[,3]==or[rpos])
  rp<- l1[indy, 1]
  rc<- l1[indy, 2]
  
  ex<- floor(nrow(dat) * (datDens[rp,rc])) #existing count of samples within the selcted quantile
  eq<- ceiling(nrow(dat) * (covDens[rp,rc])) # number of samples needed to get to equal density between data and covariates
  sn<- eq-ex #number of samples needed
  if (upSamp < sn) {sn <- upSamp} # just so we dont over allocate
  
  
  #covariate selection
  covL<- q.mat[rp, rc]
  covU<- q.mat[rp+1, rc]
  subDat<- tempD[tempD[,(base+rc)] >= covL & tempD[,(base+rc)] <= covU,] # subset the covariates that meet the standard
  
  training <- sample( nrow(subDat), sn) #random number
  subDat2<- subDat[training,]
  # subDat2 <- subDat2 %>% dplyr::select(-raw_2011_ascii,-thppm)
  
  #remove selcted samples from tempD so that repeated sampling does not occur (Is this necessary??)
  tempD<- tempD[!(tempD$cellNos %in% subDat2$cellNos), ]
  
  # Append new data to sampling dataframe
  dat<- rbind(dat,subDat2)
  
  #adjust the while params
  rpos<- rpos + 1 
  upSamp<- upSamp - sn
  print(sn)
}




# Check the sampling density with the addtional samples added

#sample data hypercube
h.mat<- matrix(0, nrow=100, ncol=7)

for (i in 1:nrow(dat)){ # the number of observations
  cntj<- 1 
  for (j in 4:ncol(dat)){ #for each column
    dd<- dat[i,j]  
    for (k in 1:100){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat[k, cntj]<- h.mat[k, cntj] + 1}
    }
    cntj<- cntj+1
  }
}

#density
h.mat[which(h.mat==0)]<- 0.000001
h.mat
datDens<-h.mat/nrow(dat)


#### check

rat<- datDens/covDens # ratio of data density and covariate density
or<- order(rat) # rank the index where the biggest discrepancy is
or

## indexes of quantiles that are not adequately sampled
l1<- which(rat < 1, arr.ind = T)
l1<- cbind(l1,which(rat < 1) )
l1
length(rat)
nrow(l1)


# What the the level of the greatest discrepancy?
indy<- which(l1[,3]==or[1])
rp<- l1[indy, 1]
rc<- l1[indy, 2]
rat[rp, rc]
h.mat[rp, rc]
cov.mat[rp, rc]
datDens[rp, rc]
covDens[rp, rc]



## The following code does not have too much to do with the algorithm

# Specify the different surveys (original and addtional)
dat$survey<- NA
str(dat)
dat[dat$cellNos==0,"survey"]<- 1
dat[dat$cellNos!=0,"survey"]<- 2

## Spatial points
coordinates(dat) <- ~x + y
str(dat)

## Coordinate reference systems
proj4string(dat) <- CRS("+init=epsg:32756")
dat@proj4string

## Write point data to shapefile
writeOGR(dat, ".", "HV_dat_shape", "ESRI Shapefile")





# The following raster has been derived by comparing the similarity between the mulivariate values of the grids and observed data
# Low numbers mean that there are not many data points showing similarity to the grid cell location. High mumber means quite a few
# obervations are similar

#raster extraction
#covariates
files<- list.files(path= getwd(), pattern = ".tif$", full.names = TRUE)
files
r1<- raster(files[1])
r1
plot(r1)

DSM_data2<- raster::extract(r1,dat, sp= 1, method = "simple")

# write table to file
write.table(as.data.frame(DSM_data2), "HELS_dat.txt", sep = ",", col.names = T)



# Doing some summary statistics between the raster grid values and the sample sites for both original and addtional data
DSM_data2$sampleNos
dat1<- DSM_data2[DSM_data2$survey ==1, ]
sum(dat1$sampleNos >= 0 & dat1$sampleNos <= 5) / nrow(dat1)
sum(dat1$sampleNos > 5 & dat1$sampleNos <= 10) / nrow(dat1)
sum(dat1$sampleNos > 10 & dat1$sampleNos <= 20) / nrow(dat1)
sum(dat1$sampleNos > 20 & dat1$sampleNos <= 40) / nrow(dat1)
sum(dat1$sampleNos > 40) / nrow(dat1)

dat2<- DSM_data2[DSM_data2$survey !=1, ]
sum(dat2$sampleNos >= 0 & dat2$sampleNos <= 5) / nrow(dat2)
sum(dat2$sampleNos > 5 & dat2$sampleNos <= 10) / nrow(dat2)
sum(dat2$sampleNos > 10 & dat2$sampleNos <= 20) / nrow(dat2)
sum(dat2$sampleNos > 20 & dat2$sampleNos <= 40) / nrow(dat2)
sum(dat2$sampleNos > 40) / nrow(dat2)

#save.image("clhs_samp.RData") #save R session



