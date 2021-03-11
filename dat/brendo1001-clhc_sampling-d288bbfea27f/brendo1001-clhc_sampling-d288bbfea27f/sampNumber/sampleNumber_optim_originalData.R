# cLHC manuscript
## Sample number optimisation
## Original data set.
## Original data set is data that has already been collected. 
## The aim of this script is to compare the collected data to the actual population in terms of covariate coverage

## Libraries
library(raster);library(rgdal);library(tripack);library(SDMTools); library(manipulate);library(clhs);library(entropy)


## Point data
p.dat<- read.table("Z:/Dropbox/2018/cLHC_samplingPAPER/gitHub/muddles/sampNumber/T1_dat.txt", header = T, sep=",") #change directory as appropriate
str(p.dat)
coordinates(p.dat)<- ~ X_REF + Y_REF 

## Raster data
cov.dat<- read.table("Z:/Dropbox/2018/cLHC_samplingPAPER/gitHub/muddles/sampNumber/T1_covs.txt", header = T, sep=",") # change directory as appropriate
str(cov.dat)
# rasterise data and stack
s1<- stack()
for (i in 3:ncol(cov.dat)){
  s1<- stack(s1, rasterFromXYZ(cov.dat[,c(1,2,i)]))}
s1

# Intersection of point data with raster data
p.dat_I<- extract(s1,p.dat, sp= 1, method = "simple")
p.dat_I<- as.data.frame(p.dat_I)
str(p.dat_I)



## comparison of population and sample distributions


## Start TEST 1
#Kullback-Leibler (KL) divergence

#Quantiles of the population
# Number of bins
nb<- 25

#quantile matrix (of the covariate data)
q.mat<- matrix(NA, nrow=(nb+1), ncol= 4)
j=1
for (i in 3:ncol(cov.dat)){ #note the index start here
  #get a quantile matrix together of the covariates
  ran1<- max(cov.dat[,i]) - min(cov.dat[,i])
  step1<- ran1/nb 
  q.mat[,j]<- seq(min(cov.dat[,i]), to = max(cov.dat[,i]), by =step1)
  j<- j+1}
q.mat

# Hypercube of population
cov.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))
for (i in 1:nrow(cov.dat)){ # the number of pixels
  cntj<- 1 
  for (j in 3:ncol(cov.dat)){ #for each column
    dd<- cov.dat[i,j]  
    for (k in 1:nb){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}
cov.mat


####Compare whole study area covariate space with the slected sample
#sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
h.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))

for (ii in 1:nrow(p.dat_I)){ # the number of observations
  cntj<- 1 
  for (jj in 4:ncol(p.dat_I)){ #for each column
    dd<- p.dat_I[ii,jj]  
    for (kk in 1:nb){  #for each quantile
      kl<- q.mat[kk, cntj] 
      ku<- q.mat[kk+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat[kk, cntj]<- h.mat[kk, cntj] + 1}
    }
    cntj<- cntj+1
  }
}
h.mat 

#Kullback-Leibler (KL) divergence
klo.1<- KL.empirical(c(cov.mat[,1]), c(h.mat[,1])) #1
klo.1
klo.2<- KL.empirical(c(cov.mat[,2]), c(h.mat[,2])) #2
klo.2
klo.3<- KL.empirical(c(cov.mat[,3]), c(h.mat[,3])) #3
klo.3
klo.4<- KL.empirical(c(cov.mat[,4]), c(h.mat[,4])) #4
klo.4
klo<- mean(c(klo.1, klo.2,klo.3,klo.4))
print(klo) # KL divergence if the existing soil sample (N=238)
#### END First Test:







#### Second Test:
##points in polygons routine

#principal component of sample
pca.s = prcomp(p.dat_I[,4:7],scale=TRUE, center=TRUE)
scores_pca1 = as.data.frame(pca.s$x)
# plot the first 2 principal components and convex hull
rand.tr<-tri.mesh(scores_pca1[,1],scores_pca1[,2])
rand.ch<-convex.hull(rand.tr, plot.it=F) #convex hull
pr_poly = cbind(x=c(rand.ch$x),y=c(rand.ch$y)) # save the convext hull vertices
plot(scores_pca1[,1], scores_pca1[,2], xlab="PCA 1", ylab="PCA 2", xlim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])),ylim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])))
lines(c(rand.ch$x,rand.ch$x[1]), c(rand.ch$y,rand.ch$y[1]),col="red",lwd=1) # draw the convex hull(domain of prediction)


# PCA prjection of population onto the pricipal components
PCA_projection<- predict(pca.s, cov.dat[,3:6]) # Project population onto sample PC
newScores = cbind(x=PCA_projection[,1],y=PCA_projection[,2]) # PC scores of projected population

#plot the polygon and all points to be checked
plot(newScores,xlab="PCA 1", ylab="PCA 2", xlim=c(min(newScores[,1:2]), max(newScores[,1:2])),ylim=c(min(newScores[,1:2]), max(newScores[,1:2])),col='black', main='convex hull of ss')
polygon(pr_poly,col='#99999990')

#create check which points fall within the polygon
specMatch = pnt.in.poly(newScores,pr_poly)
sum(specMatch$pip)/nrow(specMatch)*100 # propertion of new spectra the fall within the convex hull
points(specMatch[which(specMatch$pip==0),1:2],pch='X', col='red')
#### END Second Test:


#### Results needed from this script
# KL statistic
print(klo) # KL divergence if the existing soil sample (N=238)

# points in polygon statistic
sum(specMatch$pip)/nrow(specMatch)*100 # propertion of new spectra the fall within the convex hull






# END SCRIPT





