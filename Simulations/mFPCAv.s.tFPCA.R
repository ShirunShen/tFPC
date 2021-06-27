rm(list=ls())
######mFPCA for analyzing the simulated data#####
setwd("~/Desktop/tFPC_final/Simulations/")
library(Matrix)
library(MASS)
library(mgcv)
library(ape)
library(fields)
library(fda)   ##for basis function and penalty matrix of time


######### data generated ############
source("Simulationsettings.R")   ## ni, n, point, edges, station, etc.

set.seed(2)
n = 500  #number of time steps
ni <- sample(50:60,n, replace = TRUE)  #vector, number of points at each time step
station <- stations(n,ni)      # matrix, all data points


##PC scores
alpha1 <- ar2(0.8,0.1,sqrt(9.0),n);
alpha2 <- ar2(0.8,0.1,sqrt(4.0),n);


#alpha1 <- rnorm(n,0,sqrt(9)) #default 9
#alpha2 <- rnorm(n,0,sqrt(4)) #default 4


##the function value of each point, which is the observed value
funcvalue <- fvalue(station[1:ni[1],1],station[1:ni[1],2],1,alpha1,alpha2)
for (i in 2:n){
    funcvalue <- c(funcvalue, fvalue(station[(sum(ni[1:i-1])+1):sum(ni[1:i]),1],station[(sum(ni[1:i-1])+1):sum(ni[1:i]),2],i,alpha1,alpha2))
}

## the observed values of data points
z <- funcvalue + rnorm(length(funcvalue),0,sqrt(0.1))

load("point.dat")   #vertexes
load("edges.dat")   #edges of triangle, stored as tri vertexes,should be counter-clockwise!


###################################### tFPCA case ########################################
set.seed(1)
d = 3;   #order of Bernstein polynomials
r = 1;   #order of smoothing condition on the boundary of triangle
J = 2;   #number of principal component vectors
p = 2;   #order of AR(p)

source("../functions/tFPC/Multi_tFPCA.R")


H = smoothness(point,edges,d,r)
qrdecom = qr(t(H))
R = qr.R(qrdecom)
Q = qr.Q(qrdecom,TRUE)
nh = nrow(H)
Q2 = Q[,(qrdecom$rank+1):ncol(Q)]  ##something revised

##raw basis function value at point station[i,j]
rawbasis <- beval(point,edges,d,station[,1],station[,2])
ortho = solve(qr.R(qr(rawbasis%*%Q2)))*sqrt(nrow(rawbasis)/3) #3 is area
omat  = Q2%*%ortho
basis = rawbasis%*%omat ##exact basis we want.
#if(basis[1,1] < 0) basis = -basis
#basis = basis/sqrt((t(basis)%*%basis)[1,1])  ##important: to make basis orthonormal


Ene = energy(point,edges,d,1) ##not yet, need to consider constraints
##1 is not the same to the "r" that we defined

P = t(omat)%*%Ene%*%(omat)  ##penalty matrix corresponding to theta


####basis function for time scale
#splinebasis = create.bspline.basis(c(0,n+1), 200,names=NULL)
splinebasis = create.fourier.basis(c(0,length(ni)),nbasis=21,period=12)
basismatrix = eval.basis(1:n, splinebasis)
names(basismatrix) = NULL
####penalty matrix for time scale
Pt = fourierpen(splinebasis)

basismatrix = cbind(basismatrix,1:n/n)
Pt = cbind(Pt,0)
Pt = rbind(Pt,0)


nb = ncol(basis)
nc = ncol(basismatrix)
I_nb = diag(nb)


#penalty parameters

##### Simplex method #####
set.seed(1)
source("../functions/tFPC/tFPCA_Penaltyparameters.R")
init = c(10,1000000,20)
penalty = simplex(K_fold=5,z,point,edges,n,d,r,J,p,station,ni,basismatrix,Pt,init,"EMalgorithm")

#save.image("comparison.RData")
#penalty=list(c(5,1,20),1)

lambmus = penalty[[1]][1]
lambmut = penalty[[1]][2]
lambpc  = penalty[[1]][3]


#source("../functions/TemporalPCA/CV_tFPCA_manifoldthetab.R")
#CrossVal(K_fold=5,z,point,edges,n,d,r,J,p,station,ni,basismatrix,Pt,lambmus,lambmut,lambpc,char="EMalgorithm_tFPCA0_manifoldthetab")





#initial values
K0 = c(0.8,0.1)
sigma20 = 1
HJ0 = diag(c(1,1))
Theta0 = matrix(runif(2*nb,-1,1),ncol=2)
thetab0 = runif(nb,0,3)
thetab0 = thetab0/norm(thetab0)
thetac0 = runif(nc,0,3)
thetac0 = thetac0/norm(thetac0) #normalize thetac








Rcpp::sourceCpp("../functions/tFPC/EMalgorithm.cpp")

lambdab = 0.001 #default = 1
lambdac = 0.001
test = EMinit(z,basis,basismatrix,ni,thetab0,thetac0,P,Pt,lambdab,lambdac)
thetab0 = test$thetab
thetac0 = test$thetac

# intial values
mu0 = basis[1:ni[1],] %*% thetab0 %*% t(thetac0) %*% basismatrix[1,]
for(t in 2:n){
    temp = basis[sum(ni[1:(t-1)]) + 1:ni[t], ] %*% thetab0 %*% t(thetac0) %*% basismatrix[t,]
    mu0 = c(mu0,temp)
}
pcpart = z - mu0

f = NULL
for(t in 1:n){
    Bt = basis[indfunc(t,ni),]
    pcpart_t = pcpart[indfunc(t,ni)]
    ft = solve(t(Bt) %*% Bt +  0.001 * diag(ncol(Bt))) %*% t(Bt) %*% pcpart_t
    f = cbind(f,ft)
}

sv = svd(f,nu=J)
Theta0 = sv$u
HJ0 = diag(sv$d[1:J]^2)/length(ni)

# end intial values


system.time({parameters = EMalgorithm(z,basis,basismatrix,P,Pt,ni,HJ0,thetab0,thetac0,
    Theta0,K0,sigma20,lambmus,lambmut,lambpc)})



################ Model fitting criteria ################
thetabhat = parameters$thetab
thetachat = parameters$thetac
Thetahat  = parameters$Theta


muhat = basis[1:ni[1],] %*% thetabhat %*% t(thetachat) %*% basismatrix[1,]
for(t in 2:n){
    temp = basis[sum(ni[1:(t-1)]) + 1:ni[t], ] %*% thetabhat %*% t(thetachat) %*% basismatrix[t,]
    muhat = c(muhat, temp)
}

phi1hat = basis %*% parameters$Theta[,1]
phi2hat = basis %*% parameters$Theta[,2]
alpha1hat = parameters$alphahat[seq(1,2*n,by=2)];
alpha2hat = parameters$alphahat[seq(2,2*n,by=2)];
alpha1hat_seq = rep(alpha1hat,ni);
alpha2hat_seq = rep(alpha2hat,ni);
zhat = muhat + phi1hat * alpha1hat_seq + phi2hat * alpha2hat_seq;
residual = z - zhat;
average_residual_sq = sum(residual^2)/length(residual)

nn = 21
xn = seq(0.001,1.999,length=nn)
yn = seq(0.001,1.999,length=nn)

xx = rep(xn,nn)
yy = rep(yn, rep(nn,nn))

ind = which(xx<0.5 | xx>1.5 | yy<0.5 | yy>1.5)

station1 = cbind(xx[ind],yy[ind])

rawbasis1 <- beval(point,edges,d,station1[,1],station1[,2])

ortho1 = solve(qr.R(qr(rawbasis%*%Q2)))*sqrt(nrow(rawbasis)/3) #3 is area
omat1 = Q2 %*% ortho1
basis1 = rawbasis1 %*% omat1

#if(basis1[1,1]<0) basis1 = - basis1

realmu1 = NULL
for(t in 1:n){
    realmu1 = c(realmu1,mufunc(station1[,1],station1[,2],t))
}
realpc1 = phi1(station1[,1],station1[,2])
realpc2 = phi2(station1[,1],station1[,2])
zreal = realmu1 + kronecker(alpha1,realpc1)+ kronecker(alpha2,realpc2) ##+ rnorm(length(realmu1),0,0.04);


muhat1 = NULL
for(t in 1:n){
    muhat1 = c(muhat1, basis1 %*% thetabhat %*% t(thetachat) %*% basismatrix[t,])
}

pc1hat = basis1 %*% Thetahat[,1]
pc2hat = basis1 %*% Thetahat[,2]

zhat1 = as.vector(muhat1 + kronecker(alpha1hat, pc1hat) + kronecker(alpha2hat,pc2hat));


nnp = nrow(station1)


muhat1_tFPCA = muhat1
pc1hat_tFPCA = pc1hat
pc2hat_tFPCA = pc2hat
zhat1_tFPCA  = zhat1

#norm(pc1hat_tFPCA + realpc1,"2")

# mean function criteria
TMISE_mean_tFPCA = sum((realmu1 - muhat1)^2)/(nnp*n)
TMIAE_mean_tFPCA = sum(abs(realmu1 - muhat1))/(nnp*n)

V_pc = cbind(realpc1, realpc2)
V_pc = V_pc/norm(V_pc,"F")

Vhat_pc = -cbind(pc1hat, pc2hat)
Vhat_pc = Vhat_pc/norm(Vhat_pc,"F")

Q_v = qr.Q(qr(V_pc))#,TRUE)
Q_vhat = qr.Q(qr(Vhat_pc))#,TRUE)
rho = svd(t(Q_vhat)%*%Q_v)$d[ncol(Q_v)]
rho_V = svd(Vhat_pc %*% t(V_pc))$d[nrow(V_pc)]
angle_rho_V_tFPCA = acos(rho_V) * 180/pi

angle_rho_tFPCA = acos(rho) * 180/ pi      # angle sknewness of pc
RMSE_v_tFPCA = norm(V_pc-Vhat_pc,type="F")/sqrt(nnp)  # MSE
RMSE_Q_tFPCA = norm(Q_v - Q_vhat,type="F")/sqrt(nnp)  # MSE

# individual function criteria
TMISE_indiv_tFPCA = sum((zhat1 - zreal)^2)/(nnp*n)
TMIAE_indiv_tFPCA = sum(abs(zhat1 - zreal))/(nnp*n)





#################################### mFPCA case #######################################

set.seed(3)
d = 3;   #order of Bernstein polynomials
r = 1;   #order of smoothing condition on the boundary of triangle
J = 2;   #number of principal component vectors
source("../functions/mFPC/Multi_mFPCA.R")

#lambmu = 2     ##penalty parameter for thetamu
#lambpc = 5     ##penalty parameter for Theta


H = smoothness(point,edges,d,r)
qrdecom = qr(t(H))
R = qr.R(qrdecom)
Q = qr.Q(qrdecom,TRUE)
nh = nrow(H)
Q2 = Q[,(qrdecom$rank+1):ncol(Q)]  ##something revised
K = energy(point,edges,d,1) ##not yet, need to consider constraints
##1 is not the same to the "r" that we defined

##raw basis function value at point station[i,j]
rawbasis <- beval(point,edges,d,station[,1],station[,2])
ortho = solve(qr.R(qr(rawbasis%*%Q2)))*sqrt(nrow(rawbasis)/3) #3 is area
omat  = Q2%*%ortho
P = t(omat)%*%K%*%(omat)  ##penalty matrix corresponding to theta
basis = rawbasis%*%omat  ##exact basis we want.
#basis = -basis    ##make sure that basis[1,1] is positive


source("../functions/mFPC/EMinitial_mFPCA.R")
Rcpp::sourceCpp("../functions/mFPC/EMalg.cpp")

source("../functions/mFPC/CV_mFPCA.R")
source("../functions/mFPC/mFPCA_Penaltyparameters.R")
#penaltyPar = simplex(5,J,z,rawbasis,ni,Q2,K)


#save.image(file="comparison.RData") # to save the analysis procedure

penaltyPar = list(c(5.002511,8.374040))

lambmu = penaltyPar[[1]][1]
lambpc = penaltyPar[[1]][2]


####
set.seed(1)
EMinit = EMinitial(z,ni,basis,J)
sigma20 = 0.1 #EMinit$sigma20
D0 = diag(c(9,4)) #EMinit$D0
thetamu0 = EMinit$thetamu0
Theta0 = EMinit$Theta0


EM = EM_algorithm(z,basis,P,ni,lambmu,lambpc,sigma20,D0,thetamu0,Theta0)




#EM$sigma2
muhat = basis %*% (EM$thetamu)
#mufunc(station[,1],station[,2]) - muhat  ##estimation error of mu

#EM$alpha[1,] - alpha1
#EM$alpha[2,] - alpha2

phi1hat = basis%*%EM$BigTheta[,1]
phi2hat = basis%*%EM$BigTheta[,2]

#phi1hat - phi1(station[,1],station[,2])
#phi2hat - phi2(station[,1],station[,2])
##estimating fvalue
fhat = muhat[1:ni[1]] + basis[1:ni[1],]%*%EM$BigTheta%*%EM$alpha[,1]
for(i in 2:length(ni)){
    fhat = rbind(fhat, muhat[sum(ni[1:(i-1)])+1:ni[i]] + basis[sum(ni[1:(i-1)])+1:ni[i],]%*%EM$BigTheta%*%EM$alpha[,i])
}

#residual = funcvalue - fhat  ##between original function and estimating function
#par(mfrow=c(1,1))
#plot(residual)
#hist(residual)



########### new code August 16, 2019 ###########

nn = 21
xn = seq(0.001,1.999,length=nn)
yn = seq(0.001,1.999,length=nn)

xx = rep(xn,nn)
yy = rep(yn, rep(nn,nn))

ind = which(xx<0.5 | xx>1.5 | yy<0.5 | yy>1.5)

station1 = cbind(xx[ind],yy[ind])

rawbasis1 <- beval(point,edges,d,station1[,1],station1[,2])

ortho1 = solve(qr.R(qr(rawbasis%*%Q2)))*sqrt(nrow(rawbasis)/3) #3 is area
omat1 = Q2 %*% ortho1
basis1 = rawbasis1 %*% omat1


realmu1 = NULL
for(t in 1:n){
    realmu1 = c(realmu1,mufunc(station1[,1],station1[,2],t))
}
realpc1 = phi1(station1[,1],station1[,2])
realpc2 = phi2(station1[,1],station1[,2])
zreal = realmu1 + kronecker(alpha1,realpc1)+ kronecker(alpha2,realpc2)

muhat1 = NULL
for(t in 1:n){
    muhat1 = c(muhat1, basis1 %*% EM$thetamu)
}

pc1hat = basis1 %*% EM$BigTheta[,1]
pc2hat = basis1 %*% EM$BigTheta[,2]

zhat1 = as.vector(muhat1 + kronecker(EM$alpha[1,], pc1hat) + kronecker(EM$alpha[2,],pc2hat));

nnp = nrow(station1)
TMISE = sum((zhat1 - zreal)^2)/(nnp*n)


muhat1_mFPCA = muhat1
pc1hat_mFPCA = pc1hat
pc2hat_mFPCA = pc2hat
zhat1_mFPCA = zhat1


# mean function criteria
TMISE_mean_mFPCA = sum((realmu1 - muhat1)^2)/(nnp*n)
TMIAE_mean_mFPCA = sum(abs(realmu1 - muhat1))/(nnp*n)

V_pc = cbind(realpc1, realpc2)
Vhat_pc = cbind(pc1hat, pc2hat)

Q_v = qr.Q(qr(V_pc))#,TRUE)
Q_vhat = qr.Q(qr(Vhat_pc))#,TRUE)
rho = svd(t(Q_vhat)%*%Q_v)$d[ncol(Q_v)]

angle_rho_mFPCA = acos(rho) * 180/ pi      # angle sknewness of pc
RMSE_v_mFPCA = norm(V_pc-Vhat_pc,type="F")/sqrt(nnp)  # MSE
RMSE_Q_mFPCA = norm(Q_v - Q_vhat,type="F")/sqrt(nnp)  # MSE

# individual function criteria
TMISE_indiv_mFPCA = sum((zhat1 - zreal)^2)/(nnp*n)
TMIAE_indiv_mFPCA = sum(abs(zhat1 - zreal))/(nnp*n)


c(TMISE_mean_tFPCA, TMIAE_mean_tFPCA)
c(TMISE_mean_mFPCA, TMIAE_mean_mFPCA)

c(TMISE_indiv_tFPCA,TMIAE_indiv_tFPCA)
c(TMISE_indiv_mFPCA,TMIAE_indiv_mFPCA)

c(angle_rho_tFPCA, angle_rho_mFPCA)




save.image(file="comparison.RData") # to save the analysis procedure




setwd("~/Desktop/tFPC_final/Simulations/")
library(Matrix)
library(MASS)
library(mgcv)
library(ape)
library(fields)
library(fda)   ##for basis function and penalty matrix of time


load(file="comparison.RData")



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
fig01[ind] = -pc1hat_tFPCA
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
fig02[ind] = -pc2hat_tFPCA
image.plot(xn,yn,fig02,xlab="",ylab="",main="tFPC PC2" ,zlim=c(-1.5,1.5))
contour(xn,yn,fig02,add=TRUE)

fig02 = matrix(NA,ncol=nn,nrow=nn)
fig02[ind] = pc2hat_mFPCA
image.plot(xn,yn,fig02,xlab="",ylab="",main="mFPC PC2" ,zlim=c(-1.5,1.5))
contour(xn,yn,fig02,add=TRUE)

dev.off()





# pdf("mean_tFPCA.pdf",height=2*2,width=2*4+0.2)
# range = c(-20,40)
# par(mfrow=c(2,4),mar=c(2,4,3,1),oma=c(1,1,1,2))
# fig01 = matrix(NA,ncol=nn,nrow=nn)
# i=100
# fig01[ind] = realmu1[(i-1)*nnp+1:nnp]
# image.plot(xn,yn,fig01,xlab="",ylab="",main=paste("true mean at time",i,sep=" ") ,zlim=c(-20,40))
# contour(xn,yn,fig01,add=TRUE)
# 
# fig02 = matrix(NA,ncol=nn,nrow=nn)
# i=200
# fig02[ind] = realmu1[(i-1)*nnp+1:nnp]
# image.plot(xn,yn,fig02,xlab="",ylab="",main=paste("true mean at time",i,sep=" ") ,zlim=c(-20,40))
# contour(xn,yn,fig02,add=TRUE)
# 
# fig03 = matrix(NA,ncol=nn,nrow=nn)
# i=300
# fig03[ind] = realmu1[(i-1)*nnp+1:nnp]
# image.plot(xn,yn,fig03,xlab="",ylab="",main=paste("true mean at time",i,sep=" ") ,zlim=c(-20,40))
# contour(xn,yn,fig03,add=TRUE)
# 
# fig04 = matrix(NA,ncol=nn,nrow=nn)
# i=400
# fig04[ind] = realmu1[(i-1)*nnp+1:nnp]
# image.plot(xn,yn,fig04,xlab="",ylab="",main=paste("true mean at time",i,sep=" ") ,zlim=c(-20,40))
# contour(xn,yn,fig04,add=TRUE)
# 
# 
# fig11 = matrix(NA,ncol=nn,nrow=nn)
# i=100
# fig11[ind] = muhat1_tFPCA[(i-1)*nnp+1:nnp]
# image.plot(xn,yn,fig11,xlab="",ylab="",main=paste("tFPCA mean at time",i,sep=" ") ,zlim=c(-20,40))
# contour(xn,yn,fig11,add=TRUE)
# 
# fig12 = matrix(NA,ncol=nn,nrow=nn)
# i=200
# fig12[ind] = muhat1_tFPCA[(i-1)*nnp+1:nnp]
# image.plot(xn,yn,fig12,xlab="",ylab="",main=paste("tFPCA mean at time",i,sep=" ") ,zlim=c(-20,40))
# contour(xn,yn,fig12,add=TRUE)
# 
# fig13 = matrix(NA,ncol=nn,nrow=nn)
# i=300
# fig13[ind] = muhat1_tFPCA[(i-1)*nnp+1:nnp]
# image.plot(xn,yn,fig13,xlab="",ylab="",main=paste("tFPCA mean at time",i,sep=" ") ,zlim=c(-20,40))
# contour(xn,yn,fig13,add=TRUE)
# 
# fig14 = matrix(NA,ncol=nn,nrow=nn)
# i=400
# fig14[ind] = muhat1_tFPCA[(i-1)*nnp+1:nnp]
# image.plot(xn,yn,fig14,xlab="",ylab="",main=paste("tFPCA mean at time",i,sep=" ") ,zlim=c(-20,40))
# contour(xn,yn,fig14,add=TRUE)
# 
# dev.off()
# 


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




