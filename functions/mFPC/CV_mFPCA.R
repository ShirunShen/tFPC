CrossVal <- function(K_fold,J,z,rawbasis,ni,Q2,K,lambmu,lambpc){
  indseq = sample(1:ni[1], ni[1], replace=FALSE)
  for(i in 2:length(ni)){
    ind = sum(ni[1:(i-1)]) + sample(1:ni[i], ni[i], replace=FALSE)
    indseq = c(indseq,ind)
  }
  neglik = 0.0
  
  
  for(k in 1:K_fold){
    
    bindwidth = ceiling(ni[1]/K_fold)
    ni_negk = ni[1] - bindwidth
    if(k==K_fold){
      ni_negk = (K_fold-1)*bindwidth
    }
    ni_k = bindwidth
    if(k==K_fold){
      ni_k = ni[1] - (k-1)*bindwidth
    }
    
    indi = indseq[1:ni[1]]
    ind_negk = indi[-(((k-1)*bindwidth+1):(k*bindwidth))]
    ind_k = indi[((k-1)*bindwidth+1):(k*bindwidth)]
    if(k==K_fold){
      ind_k = indi[((k-1)*bindwidth+1):length(indi)]
    }
    
    
    for(j in 2:length(ni)){
      bindwidth = ceiling(ni[j]/K_fold)
      
      indi = indseq[sum(ni[1:(j-1)])+1:ni[j]]
      ni_negkj = ni[j]-bindwidth
      if(k==K_fold){
        ni_negkj = (K_fold-1)*bindwidth
      }
      
      ind_negk = c(ind_negk,indi[-(((k-1)*bindwidth+1):(k*bindwidth))])
      
      ind_kj = indi[(((k-1)*bindwidth+1):(k*bindwidth))]
      if(k==K_fold){
        ind_kj = indi[((k-1)*bindwidth+1):length(indi)]
      }
      ind_k = c(ind_k, ind_kj)
      ni_negk = c(ni_negk,ni_negkj)
      ni_kj = bindwidth
      if(k==K_fold){
        ni_kj = ni[j] - (k-1)*bindwidth
      }
      ni_k = c(ni_k, ni_kj)
    }
    z_negk = z[ind_negk]
    z_k = z[ind_k]
    #station_sub = station[ind_negk,]
    
    ##raw basis function value at point station[i,j]
    rawbasis_negk <- rawbasis[ind_negk,]
    
    rawbasis_k <- rawbasis[ind_k,]
    
    ortho = solve(qr.R(qr(rawbasis_negk%*%Q2)))*sqrt(nrow(rawbasis_negk)/3) #3 is area
    omat  = Q2%*%ortho
    P = t(omat)%*%K%*%(omat)  ##penalty matrix corresponding to theta
    basis_negk = rawbasis_negk%*%omat  ##exact basis we want.
    basis_k = rawbasis_k%*%omat
    if(basis_negk[1,1]<0) {basis_negk = -basis_negk}    ##make sure that basis[1,1] is positive
    if(basis_k[1,1]<0) {basis_k = -basis_k}
    
    EMinit = EMinitial(z_negk,ni_negk,basis_negk,J)
    sigma20 = EMinit$sigma20
    D0 = EMinit$D0
    thetamu0 = EMinit$thetamu0
    Theta0 = EMinit$Theta0
    EM = EM_algorithm(z_negk,basis_negk,P,ni_negk,lambmu,lambpc,sigma20,D0,thetamu0,Theta0)
    
    thetamuhat = EM$thetamu
    Thetahat = EM$BigTheta
    alphahat = EM$alpha
    sigma2hat = EM$sigma2
    Dhat = EM$D
    
    #likelihood = 
    lik_k = 0.0
    for(i in 1:length(ni_negk)){
      if(i==1){
        Bi = basis_k[1:ni_k[1],]
        zi = z_k[1:ni_k[1]]
      }else{
        Bi = basis_k[sum(ni_k[1:(i-1)])+1:ni_k[i],]
        zi = z_k[sum(ni_k[1:(i-1)])+1:ni_k[i]]
      }
      temp = Bi %*% Thetahat %*% Dhat %*% t(Thetahat) %*% t(Bi)
      temp = temp + diag(sigma2hat, ni_k[i])
      temp2 = svd(temp)
      Di = diag(temp2$d)
      Ui = temp2$u
      Vi = temp2$v
      temp3 = zi-Bi%*%thetamuhat
      lik_ki = as.numeric(log(sum(temp2$d))+t(temp3) %*% Vi %*% Di %*% t(Ui) %*% temp3)
      lik_k = lik_k + lik_ki
    }
    neglik = neglik + lik_k
  } #loop k
  
  print(c("the criterion of mFPCA is",neglik/sum(ni)))
  return(neglik/(sum(ni)))
}

