#!/usr/bin/env Rscript

## THIS FILE MUST BE EXECUTABLE
## chmod +x file
##

#make sure plotting in R works on servers with no graphical output
options(bitmapType='cairo')

##read parameters as described in the bash script
arguments <- commandArgs(trailingOnly=TRUE)
Ts= as.numeric(arguments[1])
tag= as.character(arguments[2])
fileinEpi= as.character(arguments[3])
fileinTemp= as.character(arguments[4])
fileinPrec= as.character(arguments[5])
daysmoothprec= as.numeric(arguments[6])
sides= as.numeric(arguments[7])
smoothTemp= as.numeric(arguments[8])

##define template for output file's names
outfile1<- paste("input/",tag,".4MCMC.data",sep="")
outfile2<- paste("input/",tag,".4MCMC.precDist.data",sep="")

##the model time step must be a full divisor or 1, such that days
##that have been observed in the data will match (integer) days
##in the simulations; check if that is the case here
if((1%%Ts)!=0){
	print("***********************************************")
	print("timestep must be a full divisor of 1")
	print(paste("timestep now is:",Ts,"remainder to 1:",1%%Ts))
	print("***********************************************")
	error("Abort.")
}

##open data sources / input
epiData=read.csv(paste("input/",fileinEpi,sep=""),sep=",",header=TRUE)
tempData=read.csv(paste("input/",fileinTemp,sep=""),sep=",",header=TRUE)
precData=read.csv(paste("input/",fileinPrec,sep=""),sep=",",header=TRUE)
print("<R> read data, sort out climatic variables...")

##memorize the time limits
minTimeUnit= min(tempData$day)
maxTimeUnit= max(tempData$day)

##create a time scale with the step of the model
simTimeScale= seq(minTimeUnit,maxTimeUnit,Ts) #to allow to guess time step index

##transform the temperature data to the time scale of the model by interpolation
interpolator= approxfun(tempData$day, tempData$minTemp, rule = 2)
temperature_t= interpolator(simTimeScale)
temperature_t[which(temperature_t<=0)]= 2.0 #avoid zeros and negatives on temperature
##make decisions and smooth temperature if wanted
temp_smooth= 1
if(smoothTemp) temp_smooth= daysmoothprec
nstepsFilter= (1/Ts)*(temp_smooth)
temperature_smooth_t= filter(temperature_t, rep(1/nstepsFilter, nstepsFilter), sides=sides, circular=TRUE)

##do the same for the series of observed cummulative case counts
interpolator<- approxfun(epiData$time, epiData$cases_cm, rule = 1)
ip_epiData_cases_cm<- interpolator(simTimeScale)
ip_epiData_cases_cm[which(is.na(ip_epiData_cases_cm))]<-0

##do the same for precipitation
interpolator= approxfun(precData$day, precData$prec, rule = 2)
prec_t= interpolator(simTimeScale)
prec_t[which(prec_t<=2)]=0 #avoid zeros and/or negatives on prec
nstepsFilter= (1/Ts)*(daysmoothprec)
prec_smooth_t= filter(prec_t, rep(1/nstepsFilter, nstepsFilter), sides=sides, circular=TRUE)
zTrans_prec_t= (prec_smooth_t-min(prec_smooth_t))/ (max(prec_smooth_t)-min(prec_smooth_t))

##do the same for humidity
interpolator= approxfun(precData$day, precData$um, rule = 2)
um_t= interpolator(simTimeScale)
um_t[which(um_t<=2)]=0 #avoid zeros and/or negatives on prec
nstepsFilter= (1/Ts)*(daysmoothprec)
um_smooth_t= filter(um_t, rep(1/nstepsFilter, nstepsFilter), sides=sides, circular=TRUE)
zTrans_um_t= (um_smooth_t-min(um_smooth_t))/ (max(um_smooth_t)-min(um_smooth_t))

##create a series that marks the exact days / timepoints that were observed in the data
known=rep(0,length(simTimeScale))
imatch_sim= which(simTimeScale %in% epiData$time)
imatch_dat= which(epiData$time %in% simTimeScale)
epiKnown= epiData$known
known[imatch_sim]= epiKnown[imatch_dat]

## export all series needed by the model, now interpolated to the model's time scale
output= data.frame(simTimeScale, temperature_smooth_t, ip_epiData_cases_cm, known, zTrans_prec_t, zTrans_um_t)
write.table(output, outfile1, col.names = FALSE, row.names = FALSE)

##export transformed means of precipitation and humidity
points= c(mean(zTrans_prec_t),sd(zTrans_prec_t),mean(zTrans_um_t),sd(zTrans_um_t))
output= matrix(points, ncol=2)
print(paste("meanPrec:", mean(zTrans_prec_t)))
print(paste("meanHumi:", mean(zTrans_um_t)))
write.table(output, outfile2, col.names = FALSE, row.names = FALSE)
