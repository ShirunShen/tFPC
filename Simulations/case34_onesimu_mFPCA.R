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

  #n = 500  #number of time steps
  #ni <- sample(50:60,n, replace = TRUE)  #vector, number of points at each time step
  ni <- sample(ni1:ni2,n,replace=TRUE)
  station <- stations(n,ni)      # matrix, all data points
  
  ##PC scores
  
  #alpha1 <- ar2(k1,k2,sqrt(HJ1),n);
  #alpha2 <- ar2(k1,k2,sqrt(HJ2),n);

  
  alpha1 <- rnorm(n,0,sqrt(HJ1));
  alpha2 <- rnorm(n,0,sqrt(HJ2));
  
  
  # # test
  # d = 3;   #order of Bernstein polynomials
  # r = 1;   #order of smoothing condition on the boundary of triangle
  # J = 2;   #number of principal component vectors
  # 
  # source("../functions/mFPC/Multi_mFPCA.R")
  # load("point.dat")   #vertexes
  # load("edges.dat")   #edges of triangle, stored as tri vertexes,should be counter-clockwise!
  # 
  # H = smoothness(point,edges,d,r)
  # qrdecom = qr(t(H))
  # R = qr.R(qrdecom)
  # Q = qr.Q(qrdecom,TRUE)
  # nh = nrow(H)
  # Q2 = Q[,(qrdecom$rank+1):ncol(Q)]  ##something revised
  # 
  # ##raw basis function value at point station[i,j]
  # rawbasis <- beval(point,edges,d,station[,1],station[,2])
  # ortho = solve(qr.R(qr(rawbasis%*%Q2)))*sqrt(nrow(rawbasis)/3) #3 is area
  # omat  = Q2%*%ortho
  # basis = rawbasis%*%omat  ##exact basis we want.
  # spn = ncol(basis)
  # thetab_test = c(1,1,1,rep(0,spn-3))
  # theta1_test = c(1,rep(0,spn-1))
  # theta2_test = c(0,1,rep(0,spn-2))
  # 
  # funcvalue <- basis[1:ni[1],] %*% thetab_test + basis[1:ni[1],] %*% theta1_test * alpha1[1] + basis[1:ni[1],] %*% theta2_test * alpha2[1]
  # for(i in 2:n){
  #   funcvalue <- c(funcvalue,basis[(sum(ni[1:i-1])+1):sum(ni[1:i]),] %*% thetab_test
  #                  + basis[(sum(ni[1:i-1])+1):sum(ni[1:i]),] %*% theta1_test * alpha1[i]
  #                  + basis[(sum(ni[1:i-1])+1):sum(ni[1:i]),] %*% theta2_test * alpha2[i])
  # }
  
  ##the function value of each point, which is the observed value
  funcvalue <- fvalue(station[1:ni[1],1],station[1:ni[1],2],1,alpha1,alpha2)
  for (i in 2:n){
      funcvalue <- c(funcvalue, fvalue(station[(sum(ni[1:i-1])+1):sum(ni[1:i]),1],station[(sum(ni[1:i-1])+1):sum(ni[1:i]),2],i,alpha1,alpha2))
  }
  
  ## the observed values of data points
  z <- funcvalue + rnorm(length(funcvalue),0,sqrt(var_noise))
  
  load("point.dat")   #vertexes
  load("edges.dat")   #edges of triangle, stored as tri vertexes,should be counter-clockwise!
  
  
  
  ######### parameters estimation #########
  source("../functions/mFPC/Multi_mFPCA.R")
  
  
  d = 3;   #order of Bernstein polynomials
  r = 1;   #order of smoothing condition on the boundary of triangle
  J = 2;   #number of principal component vectors
  
  #lambmu = 2     ##penalty parameter for thetamu
  #lambpc = 5     ##penalty parameter for Theta
  
  
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
  basis = rawbasis%*%omat  ##exact basis we want.
  #basis = -basis    ##make sure that basis[1,1] is positive

  
  Ene = energy(point,edges,d,1) ##not yet, need to consider constraints
  ##1 is not the same to the "r" that we defined


  P = t(omat)%*%Ene%*%(omat)  ##penalty matrix corresponding to theta

  

  
  source("../functions/mFPC/EMinitial_mFPCA.R")
  Rcpp::sourceCpp("../functions/mFPC/EMalg.cpp")
  
  
  ####tuning parameters for penalty####
  source("../functions/mFPC/CV_mFPCA.R")
  # test = c(1,100)
  # result = CrossVal(K_fold=5,J,z,rawbasis,ni,Q2,Ene,test[1],test[2])
  # result
  
  grid_mu = c(10,100,1000,10000)
  grid_pc = c(10,100,1000,10000)

  grid_mupc = expand.grid(grid_mu,grid_pc)
  result_penalty=NULL
  for(i in 1:nrow(grid_mupc)){
    result_penalty = c(result_penalty,CrossVal(K_fold=5,J,z,rawbasis,ni,Q2,Ene,grid_mupc[i,1],grid_mupc[i,2]))
  }
  init = as.numeric(unlist(grid_mupc[which.min(result_penalty),]))

  #init = c(1000,1000)
  source("../functions/mFPC/mFPCA_Penaltyparameters.R")
  penaltyPar = simplex(5,J,z,rawbasis,ni,Q2,Ene,init)

  
  #init = c(10,10)
  #penaltyPar = list(c(0,0))

  
  lambmu = penaltyPar[[1]][1]
  lambpc = penaltyPar[[1]][2]
  
  #if(lambpc > 3){
  #  lambpc = 3
  #}
  #if(lambmu < 100){
  #  lambmu = 300
  #}
  
  
  ####
  EMinit = EMinitial(z,ni,basis,J)
  sigma20 = EMinit$sigma20
  D0 = EMinit$D0
  thetamu0 = EMinit$thetamu0
  
  #thetamu0 = rep(0,length(thetamu0))
  
  Theta0 = EMinit$Theta0
  
  EM = EM_algorithm(z,basis,P,ni,lambmu,lambpc,sigma20,D0,thetamu0,Theta0)
  
  
  sigma2 = EM$sigma2
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

  
  # # test
  # realmu1 = NULL
  # for(t in 1:n){
  #   realmu1 = c(realmu1, basis1 %*% thetab_test)
  # }
  # realpc1 = basis1 %*% theta1_test
  # realpc2 = basis1 %*% theta2_test
  
   realmu1 = NULL
   for(t in 1:n){
       realmu1 = c(realmu1,mufunc(station1[,1],station1[,2],t))
   }
   realpc1 = phi1(station1[,1],station1[,2])
   realpc2 = phi2(station1[,1],station1[,2])

  
  realpcpart = kronecker(alpha1,realpc1) + kronecker(alpha2,realpc2)
  zreal = realmu1 + realpcpart
  
  muhat1 = NULL
  for(t in 1:n){
      muhat1 = c(muhat1, basis1 %*% EM$thetamu)
  }
  
  pc1hat = basis1 %*% EM$BigTheta[,1]
  pc2hat = basis1 %*% EM$BigTheta[,2]
  
  pcparthat = kronecker(EM$alpha[1,], pc1hat) + kronecker(EM$alpha[2,],pc2hat)
  zhat1 = as.vector(muhat1 + pcparthat);
  
  nnp = nrow(station1)
  TMISE = sum((zhat1 - zreal)^2)/(nnp*n)

  
  # mean function criteria
  TMISE_mean = sum((realmu1 - muhat1)^2)/(nnp*n)
  TMIAE_mean = sum(abs(realmu1 - muhat1))/(nnp*n)
  
  TMIAE_pcpart = sum(abs(realpcpart - as.vector(pcparthat)))/(nnp*n)
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
  
  angle_rho_mFPCA = acos(abs(sum(Vp * Vphat))) * 180/pi
  

  
  MAEalpha1 = min(mean(abs(alpha1 - EM$alpha[1,])),mean(abs(alpha1+EM$alpha[1,])))
  MAEalpha2 = min(mean(abs(alpha2 - EM$alpha[2,])),mean(abs(alpha2+EM$alpha[2,])))
  
  
  results = c(sigma2, EM$D[1,1], EM$D[2,2], angle_rho,  TMIAE_mean, TMIAE_indiv,TMIAE_pcpart,EM$iter,MAEalpha1, MAEalpha2, lambmu,lambpc,seed)
  
  as.matrix(results)
  
  str2 = paste("./results_final/",casetype,"T",n,"ni", ni1,ni2, "K",k1,k2,"HJ",HJ1,HJ2,"noise",var_noise,"ParamEstDetail.txt",sep="")
  write.table(t(results),file=str2,col.names=FALSE,row.names=FALSE,append = TRUE)
  
  result_list = list(param = results, 
                     trueMean = realmu1, 
                     truePC1 = realpc1, 
                     truePC2 = realpc2, 
                     trueIndiv = zreal, 
                     trueAlpha1 = alpha1, 
                     trueAlpha2 = alpha2, 
                     estMean = muhat1, 
                     estPC1 = pc1hat, 
                     estPC2 = pc2hat, 
                     estIndiv = zhat1, 
                     estAlpha1 = EM$alpha[1, ], 
                     estAlpha2 = EM$alpha[2, ], 
                     z = z, 
                     ni = ni, 
                     station = station)
  return(result_list)
}

n = 500
ni1 = 50
ni2 = 60
k1 = 0.8
k2 = 0.1
HJ1 = 1.0
HJ2 = 0.1
var_noise = 1.0

#casetype="mu1withAR_mFPCA"
casetype="mu1withiid_mFPCA"


library(doSNOW)
numCores = 36
cl = parallel::makeCluster(numCores)
registerDoSNOW(cl)
results = foreach(seed=1:100,.noexport="EMalgorithm") %dopar% onesimu(seed,n,ni1,ni2,k1,k2,HJ1,HJ2,var_noise,casetype)
stopCluster(cl)


params = NULL
for(i in 1:length(results)){
  params = rbind(params,results[[i]]$param)
}
res = rbind(apply(params,2,mean),apply(params,2,sd))
str1 = paste("./results_final/",casetype,"T",n,"ni", ni1,ni2, "K",k1,k2,"HJ",HJ1,HJ2,"noise",var_noise,"Summary.txt",sep="")
write.table(res,file=str1,col.names=FALSE,row.names=FALSE,append = TRUE)

save.image(file=paste("./results_final/",casetype,"T",n,"ni", ni1,ni2, "K",k1,k2,"HJ",HJ1,HJ2,"noise",var_noise,"results.RData",sep=""))
