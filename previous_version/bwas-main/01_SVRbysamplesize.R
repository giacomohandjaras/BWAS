#!/usr/bin/env Rscript
#

#example set up and call#####
###params
#savedir, dir to save resulting .Rdata outputs (character)
#sample size, sample size to run train data on (numeric)
#n iter, number of resampling iterations (numeric)
#y, character variable specifiying outcome variable to be predicted (must be in both train and test datasets)
#ncores, numeric variable specifying number of cores to use if paralleziation is performed in R environment
#unifeatselect, logical specifying whether univariate feature selection should be used
#nfeatures, if unifeatselect true, numeric variable specifying how many features retained (based on univariate correlation in train data)
#PCA, logical specifying if PCA should be used prior to model fitting
#propvarretain, numeric variable specifying proportion of variance (e.g., .5) to threshold number of components to retain
#tune, logical specifying whether hyperparameter tuning should occur
#tunefolds, numeric specifying number of tuning folds for nested cross-validation
#savechunksize, number of iterations that should be partitioned and saved as chunks 
#jobname, character specifying suffix for saving data

###the following parameters are passed via cluster-specific parallization in bash####
##load data##########################
batch<-as.numeric(Sys.getenv("BATCH"))
savedir<-as.character(Sys.getenv("savedir"))
samplesize<-as.numeric(Sys.getenv("samplesize"))
niter<-as.numeric(Sys.getenv("niter"))
y<-as.character(Sys.getenv("y"))
ncores<-as.numeric(Sys.getenv("ncores"))
unifeatselect<-as.logical(as.numeric(Sys.getenv("unifeatselect")))
nfeatures=as.numeric(Sys.getenv("nfeatures"))
PCA<-as.logical(as.numeric(Sys.getenv("PCA")))
propvarretain<-as.numeric(Sys.getenv("propvarretain"))
tune<-as.logical(as.numeric(Sys.getenv("tune")))
tunefolds<-as.numeric(Sys.getenv("tunefolds"))
savechunksize<-as.numeric(Sys.getenv("savechunksize"))
jobname<-as.character(Sys.getenv("jobname"))

print(ls())
print(Sys.time())

print(savedir)
print(samplesize)
print(niter)
print(y)
print(ncores)
print(unifeatselect)
print(nfeatures)
print(PCA)
print(propvarretain)
print(tune)
print(tunefolds)
print(savechunksize)


load("/home/bct16/data/traindatadf.Rdata")
load("/home/bct16/data/testdatadf.Rdata")

source("/home/bct16/scripts/ABCD_MP_SVR/0x_svrfuncs.R")
xcols<-grep("edge",names(traindf),value=TRUE)

traindf<-traindf[,!grepl("Vertex",names(traindf))]
testdf<-testdf[,!grepl("Vertex",names(testdf))]

gc()
print(ls())

####test run########

print(Sys.time())
svm_samplesizewrapper(traindf=traindf,testdf=testdf,y=y,xcols=xcols,tune=tune,tunefolds=tunefolds,outerfoldcores=ncores,tunecores=ncores,weights=TRUE,samplesizes=samplesize,niter=niter,unifeatselect=unifeatselect,nfeatures=nfeatures,PCA=PCA,propvarretain=propvarretain,savedir=savedir,savechunksize=savechunksize,validationcores=ncores,jobname=jobname)
print(Sys.time())



