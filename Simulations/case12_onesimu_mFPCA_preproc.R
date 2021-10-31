rm(list=ls())
setwd("~/Desktop/tFPC_final/Simulations/")

onesimu <- function(seed, n, ni1, ni2, k1, k2, HJ1, HJ2, var_noise, casetype){
  set.seed(seed)
  setwd("~/Desktop/tFPC_final/Simulations/")
  library(Matrix)
  library(MASS)
  library(mgcv)
  library(ape)
  library(fields)
  library(fda)
  
  source("case12_Simulationsettings.R")
  
  ni <- sample(ni1:ni2, n, replace = TRUE)
  station <- stations(n, ni)
  
  
  alpha1 <- ar2(k1, k2, sqrt(HJ1), n)
  alpha2 <- ar2(k1, k2, sqrt(HJ2), n)
  
  #alpha1 <- rnorm(n,0,sqrt(HJ1));
  #alpha2 <- rnorm(n,0,sqrt(HJ2));
  
  funcvalue <- fvalue(station[1:ni[1], 1], station[1:ni[1],2], 1, alpha1, alpha2)
  mureal <- mufunc(station[1:ni[1], 1], station[1:ni[1], 2],1)
  
  for (i in 2:n) {
    funcvalue <- c(funcvalue, fvalue(station[(sum(ni[1:i-1]) + 1):sum(ni[1:i]), 1], station[(sum(ni[1:i - 1]) + 1):sum(ni[1:i]), 2], i, alpha1, alpha2))
    mureal <- c(mureal, mufunc(station[(sum(ni[1:i - 1]) + 1):sum(ni[1:i]), 1], station[(sum(ni[1:i - 1]) + 1):sum(ni[1:i]), 2], i))
  }
  z <- funcvalue + rnorm(length(funcvalue), 0, sqrt(var_noise))
  
  
  load("point.dat")
  load("edges.dat")
  
  
  source("../functions/mFPC/Multi_mFPCA.R")
  
  d = 3
  r = 1
  J = 2
  H = smoothness(point, edges, d, r)
  qrdecom = qr(t(H))
  R = qr.R(qrdecom)
  Q = qr.Q(qrdecom, TRUE)
  nh = nrow(H)
  Q2 = Q[, (qrdecom$rank + 1):ncol(Q)]
  rawbasis <- beval(point, edges, d, station[, 1], station[,2])
  ortho = solve(qr.R(qr(rawbasis %*% Q2))) * sqrt(nrow(rawbasis)/3)
  omat = Q2 %*% ortho
  basis = rawbasis %*% omat
  
  
  Ene = energy(point, edges, d, 1)
  P = t(omat) %*% Ene %*% (omat)
  
  
  ### pre-processing
  tseq = rep(1:n, times = ni)
  splinebasis = create.fourier.basis(c(0, length(ni)), nbasis = 11, period = 12)
  basismatrix = eval.basis(tseq, splinebasis)
  names(basismatrix) = NULL
  basismatrix = cbind(basismatrix[], rep(1:n/n, times = ni))
  
  test_thetachat = solve(t(basismatrix) %*% basismatrix) %*% t(basismatrix) %*% z
  test_mu2hat = basismatrix %*% test_thetachat
  mu2hat = cbind(eval.basis(1:n, splinebasis), 1:n/n) %*% test_thetachat
  basis_new = apply(basis, 2, function(x) {return(x * test_mu2hat)})
  
  test_thetabhat = solve(t(basis_new) %*% basis_new) %*% t(basis_new) %*% z
  test_muhat = basis_new %*% test_thetabhat
  z_demean = z - test_muhat
  
  
  
  source("../functions/mFPC/EMinitial_mFPCA.R")
  Rcpp::sourceCpp("../functions/mFPC/EMalg.cpp")
  source("../functions/mFPC/CV_mFPCA.R")

  grid_mu = c(10, 100, 1000, 10000)
  grid_pc = c(10, 100, 1000, 10000)
  grid_mupc = expand.grid(grid_mu, grid_pc)
  result_penalty = NULL
  for (i in 1:nrow(grid_mupc)) {
    result_penalty = c(result_penalty, CrossVal(K_fold = 5,J, z, rawbasis, ni, Q2, Ene, grid_mupc[i, 1], grid_mupc[i, 2]))
  }

  init = as.numeric(unlist(grid_mupc[which.min(result_penalty), ]))


  source("../functions/mFPC/mFPCA_Penaltyparameters.R")

  penaltyPar = simplex(5, J, z, rawbasis, ni, Q2, Ene, init)


  lambmu = penaltyPar[[1]][1]
  lambpc = penaltyPar[[1]][2]
  
  
  
  EMinit = EMinitial(z_demean, ni, basis, J)
  sigma20 = EMinit$sigma20
  D0 = EMinit$D0
  thetamu0 = EMinit$thetamu0
  Theta0 = EMinit$Theta0
  thetamu0 = rep(0, length(thetamu0))
  
  
  EM = EM_algorithm(z_demean, basis, P, ni, lambmu, lambpc,sigma20, D0, thetamu0, Theta0)
  
  sigma2 = EM$sigma2
  muhat = basis %*% (EM$thetamu)
  phi1hat = basis %*% EM$BigTheta[, 1]
  phi2hat = basis %*% EM$BigTheta[, 2]
  fhat = muhat[1:ni[1]] + basis[1:ni[1], ] %*% EM$BigTheta %*% EM$alpha[, 1]
  
  for (i in 2:length(ni)) {
    fhat = rbind(fhat, muhat[sum(ni[1:(i - 1)]) + 1:ni[i]] + basis[sum(ni[1:(i - 1)]) + 1:ni[i], ] %*% EM$BigTheta %*% EM$alpha[, i])
  }
  
  
  
  nn = 21
  xn = seq(0.001, 1.999, length = nn)
  yn = seq(0.001, 1.999, length = nn)
  xx = rep(xn, nn)
  yy = rep(yn, rep(nn, nn))
  ind = which(xx < 0.5 | xx > 1.5 | yy < 0.5 | yy > 1.5)
  station1 = cbind(xx[ind], yy[ind])
  
  rawbasis1 <- beval(point, edges, d, station1[, 1], station1[,2])
  ortho1 = solve(qr.R(qr(rawbasis %*% Q2))) * sqrt(nrow(rawbasis)/3)
  omat1 = Q2 %*% ortho1
  basis1 = rawbasis1 %*% omat1
  
  
  realmu1 = NULL
  for (t in 1:n) {
    realmu1 = c(realmu1, mufunc(station1[, 1], station1[, 2], t))
  }
  
  realpc1 = phi1(station1[, 1], station1[, 2])
  realpc2 = phi2(station1[, 1], station1[, 2])
  realpcpart = kronecker(alpha1, realpc1) + kronecker(alpha2,realpc2)
  
  zreal = realmu1 + realpcpart
  
  
  muhat1 = NULL
  for (t in 1:n) {
    muhat1 = c(muhat1, basis1 %*% test_thetabhat * mu2hat[t])
  }
  pc1hat = basis1 %*% EM$BigTheta[, 1]
  pc2hat = basis1 %*% EM$BigTheta[, 2]
  pcparthat = kronecker(EM$alpha[1, ], pc1hat) + kronecker(EM$alpha[2,], pc2hat)
  
  
  zhat1 = as.vector(muhat1 + pcparthat)
  
  
  nnp = nrow(station1)
  
  TMISE = sum((zhat1 - zreal)^2)/(nnp * n)
  TMISE_mean = sum((realmu1 - muhat1)^2)/(nnp * n)
  TMIAE_mean = sum(abs(realmu1 - muhat1))/(nnp * n)
  TMIAE_pcpart = sum(abs(realpcpart - as.vector(pcparthat)))/(nnp * n)
  
  V_pc = cbind(realpc1, realpc2)
  Vhat_pc = cbind(pc1hat, pc2hat)
  Q_v = qr.Q(qr(V_pc))
  Q_vhat = qr.Q(qr(Vhat_pc))
  rho = svd(t(Q_vhat) %*% Q_v)$d[ncol(Q_v)]
  angle_rho = acos(rho) * 180/pi
  MSE_v = norm(V_pc - Vhat_pc, type = "F")
  MSE_Q = norm(Q_v - Q_vhat, type = "F")
  
  
  TMISE_indiv = sum((zhat1 - zreal)^2)/(nnp * n)
  TMIAE_indiv = sum(abs(zhat1 - zreal))/(nnp * n)
  
  
  Vp = c(realpc1, realpc2)
  Vp = Vp/norm(Vp, type = "2")
  Vphat = -c(pc1hat, pc2hat)
  Vphat = Vphat/norm(Vphat, type = "2")
  
  
  angle_rho_mFPCA = acos(abs(sum(Vp * Vphat))) * 180/pi
  
  
  results = c(sigma2, EM$D[1, 1], EM$D[2, 2], angle_rho, TMIAE_mean,TMIAE_indiv, TMIAE_pcpart, angle_rho_mFPCA, EM$iter, init[1], init[2], lambmu, lambpc, seed)
  
  as.matrix(results)
  
  str2 = paste("./results_final/", casetype, "T", n, "ni", ni1, ni2, "K", k1, k2, "HJ", HJ1, HJ2, "noise", var_noise, "ParamEstDetail.txt", sep = "")
  
  write.table(t(results), file = str2, col.names = FALSE, row.names = FALSE, append = TRUE)
  
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

casetype="mu1mu2withAR_mFPCA_preproc"
#casetype="mu1mu2withiid_mFPCA_preproc"

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
