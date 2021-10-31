CrossVal <- function(K_fold,basis,P,basismatrix,Pt,z,point,edges,n,d,r,J,p,ni,lambmus,lambmut,lambpc,char){
    source("../functions/tFPC/Multi_tFPCA.R")
    path0 = paste("../functions/tFPC/",char,".cpp",sep="")
    Rcpp::sourceCpp(path0)
    
    nb = ncol(P)
    nc = ncol(basismatrix)
    I_nb = diag(nb)
    
    set.seed(3)
    indseq = sample(1:ni[1], ni[1], replace=FALSE)
    for(i in 2:length(ni)){
        ind = sum(ni[1:(i-1)]) + sample(1:ni[i], ni[i], replace=FALSE)
        indseq = c(indseq,ind)
    }
    
    
  #ave_tMSE = 0.0
  sum_neglik = 0.0
  
  for(k in 1:K_fold){
    print(paste("currently run the ", k, "-th fold",sep=""))
    
    bindwidth = ceiling((ni[1])/K_fold)
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
      bindwidth = ceiling((ni[j])/K_fold)
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
    
    ##raw basis function value at point station[i,j]
    basis_negk <- basis[ind_negk,]

    
    basis_k <- basis[ind_k,]
    
    
    #print("ni_negk")
    #print(length(ni_negk))
    
    #initial values
    #set.seed(2)
    K0 = rep(0.5,p)
    sigma20 = 1.0
    HJ0 = diag(J)
    Theta0 = matrix(runif(J*nb,-1,1),ncol=J)
    thetab0 = runif(nb,-1,1)
    thetab0 = thetab0/norm(thetab0)
    thetac0 = runif(nc,-1,1)
    thetac0 = thetac0/norm(thetac0)
    
    
    #test, Oct 26, 2021
    lambdab = 0.001
    lambdac = 0.001
    test = EMinit(z_negk,basis_negk,basismatrix,ni_negk,thetab0,thetac0,P,Pt,lambdab,lambdac)
    thetab0 = test$thetab
    thetac0 = test$thetac
    
    mu0 = basis_negk[1:ni_negk[1],] %*% thetab0 %*% t(thetac0) %*% basismatrix[1,]
	for(t in 2:n){
    	temp = basis_negk[sum(ni_negk[1:(t-1)]) + 1:ni_negk[t], ] %*% thetab0 %*% t(thetac0) %*% basismatrix[t,]
    	mu0 = c(mu0,temp)
	}
	pcpart = z_negk - mu0

	f = NULL
	for(t in 1:n){
    	Bt = basis_negk[indfunc(t,ni_negk),]
    	pcpart_t = pcpart[indfunc(t,ni_negk)]
    	ft = solve(t(Bt) %*% Bt +  0.0001 * diag(ncol(Bt))) %*% t(Bt) %*% pcpart_t
    	f = cbind(f,ft)
	}

	sv = svd(f,nu=J)
	Theta0 = sv$u
	HJ0 = diag(sv$d[1:J]^2)/length(ni_negk)
    
    EM = EMalgorithm(z_negk,basis_negk,basismatrix,P,Pt,ni_negk,HJ0,thetab0,thetac0,Theta0,K0,sigma20,lambmus,lambmut,lambpc)
    #print("The number of iterations in EM algorithm")
    #print(EM$iter)
    
    theta_bc = EM$thetab %*% t(EM$thetac)
    Thetahat = EM$Theta
    alphahat = EM$alphahat
    sigma2hat = EM$sigma2
    HJhat = EM$HJ
    Khat  = EM$K
    
    #print("sigma2hat")
    #print(sigma2hat)

    
    loglik = 0.0
    logprior = 0.0
    for(t in 1:length(ni_k)){
      if(t==1){
        Bt_k = basis_k[1:ni_k[1],]
        zt_k = z_k[1:ni_k[1]]
      }else{
        Bt_k = basis_k[sum(ni_k[1:(t-1)])+1:ni_k[t],]
        zt_k = z_k[sum(ni_k[1:(t-1)])+1:ni_k[t]]
      }
      #At_k = Bt_k %*% kronecker(I_nb,t(basismatrix[t,]))
      alphahat_t = alphahat[(t-1)*J+1:J]
      temp = (zt_k - Bt_k %*% theta_bc %*% basismatrix[t,] - Bt_k %*% Thetahat %*% alphahat_t)
      loglik = loglik + t(temp)%*%temp #+ ni_k[t]*log(sigma2hat)
      
      
      #if(t>p){
      #  temp2 = 0.0
      #  for(l in 1:p){
      #    temp2 = temp2 + Khat[l] * alphahat[(t-l-1)*J+1:J]
      #  }
      #  logprior = logprior + t(alphahat_t - temp2)%*% HJhat %*% (alphahat_t - temp2)
      #}
    }
    #logprior = logprior + (n-p)*log(sum(sqrt(diag(HJhat))))
    
    print(loglik)
    sum_neglik = sum_neglik + loglik #+ logprior
    #print("summation of negative likelihood")
    #print(sum_neglik)
  } #loop k

  print(c("the Cross-validation result is",sum_neglik/sum(ni)))
  return(c(sum_neglik/sum(ni)))
}



