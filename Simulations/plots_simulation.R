rm(list=ls())
######mFPCA for analyzing the simulated data#####
setwd("~/Desktop/tFPC_final/Simulations/")
library(Matrix)
library(MASS)
library(mgcv)
library(ape)
library(fields)
library(fda)   ##for basis function and penalty matrix of time

load("./results_final/mu1mu2withAR_tFPCAT500ni5060K0.80.1HJ10.1noise1results.RData")

results_tFPCA = results
# each simulation results are stored in results_tFPCA[[i]]

load("./results_final/mu1mu2withAR_mFPCA_preprocT500ni5060K0.80.1HJ10.1noise1results.RData")

results_mFPCA = results

index = 1
zreal = results_tFPCA[[index]]$trueIndiv

zhat1_tFPCA = results_tFPCA[[index]]$estIndiv

zhat1_mFPCA = results_mFPCA[[index]]$estIndiv


# test only ##
load("mFPCA.RData")
mFPC = result

load("tFPCA.RData")
tFPC = result
realpc1 = mFPC$truePC1
realpc2 = mFPC$truePC2

pc1hat_tFPCA = tFPC$estPC1
pc2hat_tFPCA = tFPC$estPC2

pc1hat_mFPCA = mFPC$estPC1
pc2hat_mFPCA = mFPC$estPC2
##

realpc1 = results_tFPCA[[index]]$truePC1
realpc2 = results_tFPCA[[index]]$truePC2

pc1hat_tFPCA = results_tFPCA[[index]]$estPC1
pc2hat_tFPCA = results_tFPCA[[index]]$estPC2

pc1hat_mFPCA = results_mFPCA[[index]]$estPC1
pc2hat_mFPCA = results_mFPCA[[index]]$estPC2

realmu1 = results_tFPCA[[index]]$trueMean

muhat1_tFPCA = results_tFPCA[[index]]$estMean
muhat1_mFPCA = results_mFPCA[[index]]$estMean


#nn; nnp; xn,yn
nn=21
xn = seq(0.001,1.999,length=nn)
yn = seq(0.001,1.999,length=nn)

xx = rep(xn,nn)
yy = rep(yn, rep(nn,nn))

ind = which(xx<0.5 | xx>1.5 | yy<0.5 | yy>1.5)

nnp = length(ind)


###### plot ########
pdf("tFPCAversusmFPCA.pdf",height=2*3,width=2*4+0.6)
range = c(-20,40)
par(mfrow=c(3,4),mar=c(2,4,3,1),oma=c(1,1,1,2))
fig01 = matrix(NA,ncol=nn,nrow=nn)
i=100
fig01[ind] = zreal[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig01,xlab="",ylab="",main=paste("true function at ",i,sep=" ") ,zlim=range)
contour(xn,yn,fig01,add=TRUE)

fig02 = matrix(NA,ncol=nn,nrow=nn)
i=200
fig02[ind] = zreal[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig02,xlab="",ylab="",main=paste("true function at ",i,sep=" ") ,zlim=range)
contour(xn,yn,fig02,add=TRUE)

fig03 = matrix(NA,ncol=nn,nrow=nn)
i=300
fig03[ind] = zreal[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig03,xlab="",ylab="",main=paste("true function at ",i,sep=" ") ,zlim=range)
contour(xn,yn,fig03,add=TRUE)

fig04 = matrix(NA,ncol=nn,nrow=nn)
i=400
fig04[ind] = zreal[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig04,xlab="",ylab="",main=paste("true function at ",i,sep=" ") ,zlim=range)
contour(xn,yn,fig04,add=TRUE)


fig11 = matrix(NA,ncol=nn,nrow=nn)
i=100
fig11[ind] = zhat1_tFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig11,xlab="",ylab="",main=paste("tFPC prediction at",i,sep=" ") ,zlim=range)
contour(xn,yn,fig11,add=TRUE)

fig12 = matrix(NA,ncol=nn,nrow=nn)
i=200
fig12[ind] = zhat1_tFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig12,xlab="",ylab="",main=paste("tFPC prediction at",i,sep=" ") ,zlim=range)
contour(xn,yn,fig12,add=TRUE)

fig13 = matrix(NA,ncol=nn,nrow=nn)
i=300
fig13[ind] = zhat1_tFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig13,xlab="",ylab="",main=paste("tFPC prediction at",i,sep=" ") ,zlim=range)
contour(xn,yn,fig13,add=TRUE)

fig14 = matrix(NA,ncol=nn,nrow=nn)
i=400
fig14[ind] = zhat1_tFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig14,xlab="",ylab="",main=paste("tFPC prediction at",i,sep=" ") ,zlim=range)
contour(xn,yn,fig14,add=TRUE)



fig11 = matrix(NA,ncol=nn,nrow=nn)
i=100
fig11[ind] = zhat1_mFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig11,xlab="",ylab="",main=paste("mFPC prediction at",i,sep=" ") ,zlim=range)
contour(xn,yn,fig11,add=TRUE)

fig12 = matrix(NA,ncol=nn,nrow=nn)
i=200
fig12[ind] = zhat1_mFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig12,xlab="",ylab="",main=paste("mFPC prediction at",i,sep=" ") ,zlim=range)
contour(xn,yn,fig12,add=TRUE)

fig13 = matrix(NA,ncol=nn,nrow=nn)
i=300
fig13[ind] = zhat1_mFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig13,xlab="",ylab="",main=paste("mFPC prediction at",i,sep=" ") ,zlim=range)
contour(xn,yn,fig13,add=TRUE)

fig14 = matrix(NA,ncol=nn,nrow=nn)
i=400
fig14[ind] = zhat1_mFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig14,xlab="",ylab="",main=paste("mFPC prediction at",i,sep=" ") ,zlim=range)
contour(xn,yn,fig14,add=TRUE)

dev.off()




pdf("PC_tFPCAversusmFPCA.pdf",height=2*2,width=2*3+0.6)
par(mfrow=c(2,3),mar=c(2,4,3,1),oma=c(1,1,1,2))
fig01 = matrix(NA,ncol=nn,nrow=nn)
fig01[ind] = realpc1
image.plot(xn,yn,fig01,xlab="",ylab="",main="true PC1" ,zlim=c(-1.5,1.5))
contour(xn,yn,fig01,add=TRUE)


fig01 = matrix(NA,ncol=nn,nrow=nn)
fig01[ind] = pc1hat_tFPCA
image.plot(xn,yn,fig01,xlab="",ylab="",main="tFPC PC1" ,zlim=c(-1.5,1.5))
contour(xn,yn,fig01,add=TRUE)

fig01 = matrix(NA,ncol=nn,nrow=nn)
fig01[ind] = pc1hat_mFPCA
image.plot(xn,yn,fig01,xlab="",ylab="",main="mFPC PC1" ,zlim=c(-1.5,1.5))
contour(xn,yn,fig01,add=TRUE)


fig02 = matrix(NA,ncol=nn,nrow=nn)
fig02[ind] = realpc2
image.plot(xn,yn,fig02,xlab="",ylab="",main="true PC2" ,zlim=c(-1.5,1.5))
contour(xn,yn,fig02,add=TRUE)


fig02 = matrix(NA,ncol=nn,nrow=nn)
fig02[ind] = pc2hat_tFPCA
image.plot(xn,yn,fig02,xlab="",ylab="",main="tFPC PC2" ,zlim=c(-1.5,1.5))
contour(xn,yn,fig02,add=TRUE)

fig02 = matrix(NA,ncol=nn,nrow=nn)
fig02[ind] = pc2hat_mFPCA
image.plot(xn,yn,fig02,xlab="",ylab="",main="mFPC PC2" ,zlim=c(-1.5,1.5))
contour(xn,yn,fig02,add=TRUE)

dev.off()


pdf("mean_tFPCAversusmFPCA.pdf",height=2*3,width=2*4+0.6)
range = c(-20,40)
par(mfrow=c(3,4),mar=c(2,4,3,1),oma=c(0.5,0,0,1.5))
fig01 = matrix(NA,ncol=nn,nrow=nn)
i=100
fig01[ind] = realmu1[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig01,xlab="",ylab="",main=paste("true mean at time",i,sep=" ") ,zlim=c(-20,40))
contour(xn,yn,fig01,add=TRUE)

fig02 = matrix(NA,ncol=nn,nrow=nn)
i=200
fig02[ind] = realmu1[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig02,xlab="",ylab="",main=paste("true mean at time",i,sep=" ") ,zlim=c(-20,40))
contour(xn,yn,fig02,add=TRUE)

fig03 = matrix(NA,ncol=nn,nrow=nn)
i=300
fig03[ind] = realmu1[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig03,xlab="",ylab="",main=paste("true mean at time",i,sep=" ") ,zlim=c(-20,40))
contour(xn,yn,fig03,add=TRUE)

fig04 = matrix(NA,ncol=nn,nrow=nn)
i=400
fig04[ind] = realmu1[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig04,xlab="",ylab="",main=paste("true mean at time",i,sep=" ") ,zlim=c(-20,40))
contour(xn,yn,fig04,add=TRUE)


fig11 = matrix(NA,ncol=nn,nrow=nn)
i=100
fig11[ind] = muhat1_tFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig11,xlab="",ylab="",main=paste("tFPC mean at time",i,sep=" ") ,zlim=c(-20,40))
contour(xn,yn,fig11,add=TRUE)

fig12 = matrix(NA,ncol=nn,nrow=nn)
i=200
fig12[ind] = muhat1_tFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig12,xlab="",ylab="",main=paste("tFPC mean at time",i,sep=" ") ,zlim=c(-20,40))
contour(xn,yn,fig12,add=TRUE)

fig13 = matrix(NA,ncol=nn,nrow=nn)
i=300
fig13[ind] = muhat1_tFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig13,xlab="",ylab="",main=paste("tFPC mean at time",i,sep=" ") ,zlim=c(-20,40))
contour(xn,yn,fig13,add=TRUE)

fig14 = matrix(NA,ncol=nn,nrow=nn)
i=400
fig14[ind] = muhat1_tFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig14,xlab="",ylab="",main=paste("tFPC mean at time",i,sep=" ") ,zlim=c(-20,40))
contour(xn,yn,fig14,add=TRUE)



fig11 = matrix(NA,ncol=nn,nrow=nn)
i=100
fig11[ind] = muhat1_mFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig11,xlab="",ylab="",main=paste("mFPC mean at time",i,sep=" ") ,zlim=range)
contour(xn,yn,fig11,add=TRUE)

fig12 = matrix(NA,ncol=nn,nrow=nn)
i=200
fig12[ind] = muhat1_mFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig12,xlab="",ylab="",main=paste("mFPC mean at time",i,sep=" ") ,zlim=range)
contour(xn,yn,fig12,add=TRUE)

fig13 = matrix(NA,ncol=nn,nrow=nn)
i=300
fig13[ind] = muhat1_mFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig13,xlab="",ylab="",main=paste("mFPC mean at time",i,sep=" ") ,zlim=range)
contour(xn,yn,fig13,add=TRUE)

fig14 = matrix(NA,ncol=nn,nrow=nn)
i=400
fig14[ind] = muhat1_mFPCA[(i-1)*nnp+1:nnp]
image.plot(xn,yn,fig14,xlab="",ylab="",main=paste("mFPC mean at time",i,sep=" ") ,zlim=range)
contour(xn,yn,fig14,add=TRUE)
dev.off()




