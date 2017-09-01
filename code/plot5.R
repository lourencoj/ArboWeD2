#!/usr/bin/env Rscript

## THIS FILE MUST BE EXECUTABLE
## chmod +x file
##

#make sure plotting in R works on servers with no graphical output
options(bitmapType='cairo')

##read some parameters as defined in the bash script
arguments <- commandArgs(trailingOnly=TRUE)
eta= as.numeric(arguments[1])
tag= as.character(arguments[2])
fileinEpi= as.character(arguments[3])
tRho= as.numeric(arguments[4])
NH= as.numeric(arguments[5])

require("scales")


print("reading input for MCMC")
data=read.table(paste("input/",tag,".4MCMC.data",sep=""), header=FALSE)
colnames(data)=c("time","temp","cases_cm")

##decide on model step size
stepsize=data$time[2]-data$time[1]

print("reading original epi data")
epi=read.table(paste("input/",fileinEpi,sep=""), header=TRUE, sep=",")

print("reading accepted parameters")
res= read.table(paste("output/",tag,".accepted_parameters.data",sep=""), header=TRUE)
res.burnedin= res[(nrow(res)*(1-eta)):nrow(res),]
write.table(res.burnedin, paste("output/",tag,".accepted_parameters_burnedIn.csv",sep=""), col.names=TRUE, row.names=FALSE, sep=' ') #save the parameters after burn in
res.burnedin.mean= matrix(colMeans(res.burnedin),ncol=ncol(res.burnedin))
colnames(res.burnedin.mean)= colnames(res.burnedin)
write.table(res.burnedin.mean, paste("output/",tag,".accepted_parameters_burnedIn_mean.csv",sep=""), col.names=TRUE, row.names=FALSE, sep=' ') #save the parameters after burn in

print("reading inc_t")
inc_t= as.numeric(unlist(read.table(paste("output/",tag,".inc_t.data",sep=""), header=FALSE, col.names=FALSE)))

print("reading RH_t")
RH_t= as.numeric(unlist(read.table(paste("output/",tag,".RH_t.data",sep=""), header=FALSE, col.names=FALSE)))

print("reading v_t")
v_t= as.numeric(unlist(read.table(paste("output/",tag,".V_t.data",sep=""), header=FALSE, col.names=FALSE)))

print("reading A_t")
A_t= as.numeric(unlist(read.table(paste("output/",tag,".A_t.data",sep=""), header=FALSE, col.names=FALSE)))

print("reading r0_t")
r0_t= as.numeric(unlist(read.table(paste("output/",tag,".r0_t.data",sep=""), header=FALSE, col.names=FALSE)))

print("reading re_t")
re_t= as.numeric(unlist(read.table(paste("output/",tag,".re_t.data",sep=""), header=FALSE, col.names=FALSE)))

##read climatic series
inD= read.table(paste("input/",tag,".4MCMC.data",sep=""))
rainK_t= as.numeric(inD[,5])
humi_t= as.numeric(inD[,6])

#########################################################

filename<- paste0('output/',tag,'.weather.vs.vector.R0Re.png')
png(filename,w=1600,h=450,bg='white')
	layout(matrix(seq(1,3), ncol = 3, byrow=T))
	par(mar=c(4, 4, 2, 4), cex=1)

	Y<- as.numeric(A_t)
	yrange<-c(0,max(Y)*1.25)
	plot(data$time,Y, t='l', ylim=yrange, lwd=2, col="purple", main="weather vs A", ylab="A", xlab="days")
	par(new=TRUE)
		Y= as.numeric(rainK_t)
		Z= as.numeric(humi_t)
		plot(data$time, Y,type="l",col="orange",xaxt="n",yaxt="n",xlab="",ylab="")
		lines(data$time, Z, t='l', col="black")
		lines(data$time,data$temp/max(data$temp),col="grey", lwd=1.5)
		axis(4)
		mtext("prec / humi / temp",side=4,line=3)
		box(lwd=2)

	Y<- as.numeric(v_t)
	yrange<-c(0,max(Y)*1.25)
	plot(data$time,Y, t='l', ylim=yrange, lwd=2, col="blue", main="weather vs V", ylab="V", xlab="days")
	par(new=TRUE)
		Y= as.numeric(rainK_t)
		Z= as.numeric(humi_t)
		plot(data$time, Y,type="l",col="orange",xaxt="n",yaxt="n",xlab="",ylab="")
		lines(data$time, Z, t='l', col="black")
		lines(data$time,data$temp/max(data$temp),col="grey", lwd=1.5)
		axis(4)
		mtext("prec / humi / temp",side=4,line=3)
		box(lwd=2)

		Y<- as.numeric(r0_t)
		yrange<-c(0,max(Y)*1.25)
		plot(data$time,Y, t='l', ylim=yrange, lwd=2, col="red", main="weather vs R0 Re", ylab="R0 Re", xlab="days")
		lines(data$time,re_t,t='l',col='cyan')
		par(new=TRUE)
			Y= as.numeric(rainK_t)
			Z= as.numeric(humi_t)
			plot(data$time, Y,type="l",col="orange",xaxt="n",yaxt="n",xlab="",ylab="")
			lines(data$time, Z, t='l', col="black")
			lines(data$time,data$temp/max(data$temp),col="grey", lwd=1.5)
			axis(4)
			mtext("prec / humi / temp",side=4,line=3)
			box(lwd=2)

a=dev.off()



filename<- paste("output/",tag,".epidemic.fit.png",sep="")
png(filename, height=400, width=800, bg="white" )

	layout(matrix(seq(1,2), ncol = 2, byrow=T))
	par(mar=c(4, 4, 2, 1.15), cex=1)

	Y<- as.numeric(RH_t)
	yrange<-c(0,max(epi$cases_cm)*1.2)
	plot(data$time,Y, t='n', xaxs="i", yaxs="i", ylim=yrange, lwd=2, col="gray30", main="data vs zeta*fit", ylab="cummulative cases", xlab="days")
	lines(epi$time, epi$cases_cm, lwd=3)
	lines(data$time,Y,col="blue", lwd=1.5)
	legend('topleft',legend=c('data','model'),col=c('black','blue'),lwd=2)
	box(lwd=2)

	Y<- as.numeric(inc_t)
	yrange<-c(0,max(c(epi$cases,data$temp))*1.2)
	xrange<-c(-50,100)
	plot(data$time,Y, t='n', xaxs="i", yaxs="i", ylim=yrange, lwd=2, main="estimated real epidemic", ylab="cases", xlab="days")
	lines(data$time,Y,col="blue", lwd=1.5)
	lines(data$time,data$temp,col="grey", lwd=1.5)
	Y1<- as.numeric(v_t)
	lines(data$time,Y1, col="limegreen", t='l', lwd=2)
	Y1<- as.numeric(r0_t)
	lines(data$time,Y1, col="red", t='l', lwd=2)
	abline(h=1, col="red", lty=2)
	legend('topleft',legend=c('model epidemic','temperature','R0','V'),col=c('blue','grey','red','green'),lwd=2)
	box(lwd=2)


a=dev.off()



filename<- paste("output/",tag,".mcmc.png",sep="")
png(filename, height=600, width=900, bg="white" )

	layout(matrix(seq(1,6), ncol = 3, byrow=T))
	par(mar=c(4, 4, 1.5, 1.15), cex=1)

	ref<- round(nrow(res)*(1-eta))

	Y<- res$T0
	plot(Y, t='l', main="T0", xlab="MCMC state", ylab="T0 accepted", xaxs='i', yaxs='i', lwd=3)
	abline(v=ref)
	box(lwd=2)

	Y<- res$K
	plot(Y, t='l', main="K", xlab="MCMC state", ylab="K accepted", xaxs='i', yaxs='i', lwd=3)
	abline(v=ref)
	box(lwd=2)


	Y1<- res$ecoFactorV
	Y2<- res$eco2FactorV
	Y3<- res$eipFactorV
	yylim<- range(Y1,Y2,Y3)
	plot(Y1, t='l', main="factors", xlab="MCMC state", ylab="factors accepted", xaxs='i', yaxs='i', lwd=3, ylim=yylim)
	lines(Y2, col='brown', lwd=3)
	lines(Y3, col='blue', lwd=3)
	abline(v=ref)
	legend('bottomright',legend=c('alpha','rho','eta'),col=c('black','brown','blue'),lwd=2)
	box(lwd=2)


	Y1<- res$eipH
	Y2<- res$infpH
	yylim<- range(Y1,Y2)
	plot(Y1, t='l', main="human", xlab="MCMC state", ylab="human parms accepted", xaxs='i', yaxs='i', lwd=3, ylim=yylim)
	lines(Y2, col='brown', lwd=3)
	abline(v=ref)
	legend('bottomright',legend=c('eipH','infpH'),col=c('black','brown'),lwd=2)
	box(lwd=2)


	Y<- res$zeta
	plot(Y, t='l', main="zeta", xlab="MCMC state", ylab="zeta accepted", xaxs='i', yaxs='i', lwd=3)
	abline(v=ref)
	box(lwd=2)

	Y<- res$rho
	plot(Y, t='l', main="rho", xlab="MCMC state", ylab="rho accepted", xaxs='i', yaxs='i', lwd=3)
	abline(v=ref)
	box(lwd=2)

a=dev.off()


filename<- paste("output/",tag,".aprob.png",sep="")
png(filename, height=350, width=900, bg="white" )

	layout(matrix(seq(1,2), ncol = 2, byrow=T))
	par(mar=c(4, 4, 1.5, 1.15), cex=1)

	accepted=which(res$accept==1)
	aProbAccepted= (res$aProb)[accepted]
	plot(aProbAccepted, t='l', xaxs="i", ylim=c(0,1), main='probability of acceptance')
	box(lwd=2)

	f20 <- rep(1/100, 100)
	y_lag <- filter(res$accept, f20, sides=1)
	plot(y_lag, main="smoothed acceptance rate", ylim=c(0,1))

a=dev.off()
