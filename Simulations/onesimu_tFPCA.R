# June 16, change the qr decomposition condition, complete = FALSE

rm(list=ls())
setwd("~/Desktop/tFPC_final/Simulations/")
onesimu <- function(seed,n,ni1,ni2,k1,k2,HJ1,HJ2,var_noise,casetype){
  set.seed(seed)
  ########tFPCA for analyzing the simulated data#######
  setwd("~/Desktop/tFPC_final/Simulations/")
  library(Matrix)
  library(MASS)
  library(mgcv)
  library(ape)
  library(fields)
  library(fda)   ##for basis function and penalty matrix of time
  
  source("Simulationsettings.R")  ##data
  source("../functions/tFPC/Multi_tFPCA.R")
  
  #n = 500  #number of time steps
  #ni <- sample(50:60,n, replace = TRUE)  #vector, number of points at each time step
  ni <- sample(ni1:ni2,n,replace=TRUE)
  station <- stations(n,ni)      # matrix, all data points
  
  ##PC scores
  
  alpha1 <- ar2(k1,k2,sqrt(HJ1),n);
  alpha2 <- ar2(k1,k2,sqrt(HJ2),n);
  
  #alpha1 <- rnorm(n,0,sqrt(HJ1));
  #alpha2 <- rnorm(n,0,sqrt(HJ2));
  
  ##the function value of each point, which is the observed value
  funcvalue <- fvalue(station[1:ni[1],1],station[1:ni[1],2],1,alpha1,alpha2)
  for (i in 2:n){
      funcvalue <- c(funcvalue, fvalue(station[(sum(ni[1:i-1])+1):sum(ni[1:i]),1],station[(sum(ni[1:i-1])+1):sum(ni[1:i]),2],i,alpha1,alpha2))
  }
  
  ## the observed values of data points
  z <- funcvalue + rnorm(length(funcvalue),0,sqrt(var_noise))
  
  load("point.dat")   #vertexes
  load("edges.dat")   #edges of triangle, stored as tri vertexes,should be counter-clockwise!
  
  ##### parameter estimation #####
  
  d = 3;   #order of Bernstein polynomials
  r = 1;   #order of smoothing condition on the boundary of triangle
  J = 2;   #number of principal component vectors
  p = 2;   #order of AR(p)
  
  H = smoothness(point,edges,d,r)
  qrdecom = qr(t(H))
  R = qr.R(qrdecom)
  Q = qr.Q(qrdecom,TRUE)
  nh = nrow(H)
  Q2 = Q[,(qrdecom$rank+1):ncol(Q)]  ##something revised
  
  ####basis function for time scale
  #splinebasis = create.bspline.basis(c(0,n+1), 200,names=NULL)
  splinebasis = create.fourier.basis(c(0,length(ni)),nbasis=21,period=12) #nbasis=3 is an issue. June 10, 2020
  basismatrix = eval.basis(1:n, splinebasis)
  names(basismatrix) = NULL
  ####penalty matrix for time scale
  Pt = fourierpen(splinebasis)
  
  
  # modified, Aug 27,2020
  basismatrix = cbind(basismatrix,1:n/n)
  Pt = cbind(Pt,0)
  Pt = rbind(Pt,0)
  
  
  ##raw basis function value at point station[i,j]
  rawbasis <- beval(point,edges,d,station[,1],station[,2])
  ortho = solve(qr.R(qr(rawbasis%*%Q2)))*sqrt(nrow(rawbasis)/3) #3 is area
  omat  = Q2%*%ortho
  basis = rawbasis%*%omat ##exact basis we want.
  #basis = basis/sqrt((t(basis)%*%basis)[1,1])  ##important: to make basis orthonormal
  
  
  Ene = energy(point,edges,d,1) ##not yet, need to consider constraints
  ##1 is not the same to the "r" that we defined
  
  P = t(omat)%*%Ene%*%(omat)  ##penalty matrix corresponding to theta
  nb = ncol(basis)
  nc = ncol(basismatrix)
  I_nb = diag(nb)
  
  
  #initial values
  K0 = c(0.8,0.1)
  sigma20 = 1
  HJ0 = diag(2)
  Theta0 = matrix(runif(2*nb,-1,1),ncol=2)
  thetab0 = runif(nb,-1,1)
  thetab0 = thetab0/norm(thetab0)
  thetac0 = runif(nc,-1,1)
  thetac0 = thetac0/norm(thetac0) #normalize thetac
  
  #penalty parameters
  
  ##### without penalty #####
  
  #penalty = list(c(0.0,0.0,0.0),1)
  
  
  
  ##### Simplex method #####
  
  
  init = c(10,1,20) # default (10,1,20)
  source("../functions/tFPC/tFPCA_Penaltyparameters.R")
  penalty = simplex(K_fold=5,z,point,edges,n,d,r,J,p,station,ni,basismatrix,Pt,init,"EMalgorithm")

# penalty = list(c(8.2165658,4.3437858,7.53277),1) #seed =3
# penalty = list(c(0.8767,0,28.8654167),1) #seed = 1
# penalty = list(c(1,10,20),1)
# penalty = list(c(5.888278,8.2343645,8.095942),1)

  lambmus = penalty[[1]][1]
  lambmut = penalty[[1]][2]
  lambpc  = penalty[[1]][3]
  iters   = penalty[[2]]
  
  Rcpp::sourceCpp("../functions/tFPC/EMalgorithm.cpp")
  
  # update the initial values of thetab0 and thetac0
  lambdab = 0.001 #default = 1
  lambdac = 10000
  test = EMinit(z,basis,basismatrix,ni,thetab0,thetac0,P,Pt,lambdab,lambdac)
  thetab0 = test$thetab
  thetac0 = test$thetac
  
  # initial values
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
    ft = solve(t(Bt) %*% Bt + 0.001 *diag(ncol(Bt))) %*% t(Bt) %*% pcpart_t
    f = cbind(f,ft)
  }
  
  sv = svd(f,nu=J)
  Theta0 = sv$u
  HJ0 = diag(sv$d[1:J]^2)/length(ni)
  # end initial values
  
  parameters = EMalgorithm(z,basis,basismatrix,P,Pt,ni,HJ0,thetab0,thetac0,
      Theta0,K0,sigma20,lambmus,lambmut,lambpc)
  
  k1hat = parameters$K[1]
  k2hat = parameters$K[2]
  sigma_noise = parameters$sigma2
  sigma_eta1 = parameters$HJ[1,1]
  sigma_eta2 = parameters$HJ[2,2]
  
  
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
  zreal = realmu1 + kronecker(alpha1,realpc1)+ kronecker(alpha2,realpc2)
  
  muhat1 = NULL
  for(t in 1:n){
      muhat1 = c(muhat1, basis1 %*% thetabhat %*% t(thetachat) %*% basismatrix[t,])
  }
  
  pc1hat = basis1 %*% Thetahat[,1]
  pc2hat = basis1 %*% Thetahat[,2]
  
  zhat1 = as.vector(muhat1 + kronecker(alpha1hat, pc1hat) + kronecker(alpha2hat,pc2hat));
  
  nnp = nrow(station1)

  
  ######## model fitting above ###########
  
  
  # mean function criteria
  TMISE_mean = sum((realmu1 - muhat1)^2)/(nnp*n)
  TMIAE_mean = sum(abs(realmu1 - muhat1))/(nnp*n)
  
  
  # principal component criteria: principal angle, MSE of original PC matrix, svd pc orthonormal matrix
  V_pc = cbind(realpc1, realpc2)
  Vhat_pc = cbind(pc1hat, pc2hat)
  
  Q_v = qr.Q(qr(V_pc))#,TRUE)
  Q_vhat = qr.Q(qr(Vhat_pc))#,TRUE)
  rho = svd(t(Q_vhat)%*%Q_v)$d[ncol(Q_v)]
  
  angle_rho = acos(rho) * 180/ pi      # angle sknewness of pc
  MSE_v = norm(V_pc-Vhat_pc,type="F")  # MSE
  MSE_Q = norm(Q_v - Q_vhat,type="F")  # MSE
  
  # individual function criteria
  TMISE_indiv = sum((zhat1 - zreal)^2)/(nnp*n)
  TMIAE_indiv = sum(abs(zhat1 - zreal))/(nnp*n)
  
  
  Vp = c(realpc1, realpc2)
  Vp = Vp/norm(Vp,type="2")
  Vphat = c(pc1hat,pc2hat)
  Vphat = Vphat/norm(Vphat,type="2")
  
  angle_rho_tFPCA = acos(abs(sum(Vp * Vphat))) * 180/pi
  

  
  results = c(k1hat,k2hat,sigma_noise,sigma_eta1,sigma_eta2,TMISE_mean,TMIAE_mean, angle_rho, MSE_v, MSE_Q, TMISE_indiv, TMIAE_indiv ,iters, lambmus, lambmut, lambpc,angle_rho_tFPCA,as.numeric(seed))
  
  str2 = paste("./results_tFPCA/1",casetype,"T",n,"ni", ni1,ni2, "K",k1,k2,"HJ",HJ1,HJ2,"noise",var_noise,"ParamEstDetail.txt",sep="")
  write.table(t(results),file=str2,col.names=FALSE,row.names=FALSE,append = TRUE)
  
  #print("tFPCA")
  #print(results)
  
  return(results)
}

n = 500
ni1 = 50
ni2 = 60
k1 = 0.8
k2 = 0.1
HJ1 = 49.0
HJ2 = 25.0
var_noise = 1.0


#casetype = "mu1mu2withAR"
casetype="mu1withAR_penalty"
#casetype="mu1mu2withiid"
#casetype="mu1withiid"

library(doSNOW)
numCores = 32
cl = parallel::makeCluster(numCores)
registerDoSNOW(cl)
results = foreach(seed=1:100,.combine=rbind,.noexport="EMalgorithm") %dopar% onesimu(seed,n,ni1,ni2,k1,k2,HJ1,HJ2,var_noise,casetype)
stopCluster(cl)


#casetype="magnitude10_HJ51"


res = rbind(apply(results,2,mean),apply(results,2,sd))
#res = summary(results)
str1 = paste("./results_tFPCA/1",casetype,"T",n,"ni", ni1,ni2, "K",k1,k2,"HJ",HJ1,HJ2,"noise",var_noise,"Summary.txt",sep="")
write.table(res,file=str1,col.names=FALSE,row.names=FALSE,append = TRUE)
