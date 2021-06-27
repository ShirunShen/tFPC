AIC <- function(p,K_fold,basis,P,basismatrix,Pt,z,point,edges,n,d,r,J,ni){


    print(p)
    
    source("../functions/tFPC/tFPCA_simplex_realdata.R")
    
    init = c(10000,0.0001,10000)
    penalty = simplex(K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,p,ni,init,"EMalgorithm_tFPCA0_manifoldthetab")
    
    #penalty = list(c(10000,0.0001,10000),1)
    lambmus = penalty[[1]][1]
    lambmut = penalty[[1]][2]
    lambpc  = penalty[[1]][3]
    
    set.seed(1)
    K0 = rep(0.5,p)
    sigma20 = 1.0
    HJ0 = diag(J)
    Theta0 = matrix(runif(J*nb,-1,1),ncol=J)
    thetab0 = runif(nb,-1,1)
    thetac0 = runif(nc,-1,1)
    thetac0 = thetac0/norm(thetac0)
    
    print("EM start")
    
    Rcpp::sourceCpp("../functions/tFPC/EMalgorithm.cpp")
    parameters = EMalgorithm(z,basis,basismatrix,P,Pt,ni,HJ0,thetab0,thetac0,
    Theta0,K0,sigma20,lambmus,lambmut,lambpc)
    
    print("EM end")
    
    thetabhat = parameters$thetab
    thetachat = parameters$thetac
    Thetahat = parameters$Theta
    alphahat = parameters$alphahat
    sigma2hat = parameters$sigma2
    HJhat = parameters$HJ
    Khat = parameters$K
    
    negloglik = 0.0
    for(t in 1:length(ni)){
        if(t==1){
            Bt = basis[1:ni[1],]
            zt = z[1:ni[1]]
        }else{
            Bt = basis[sum(ni[1:(t-1)])+1:ni[t],]
            zt = basis[sum(ni[1:(t-1)]) + 1:ni[t]]
        }
        alphahat_t = alphahat[(t-1)*J+1:J]
        temp = (zt - Bt %*% thetabhat %*% t(thetachat) %*% basismatrix[t,] - Bt %*% Thetahat %*% alphahat_t)
        negloglik = negloglik + t(temp) %*% temp / sigma2hat
    }
    
    print("AIC calculation")
    
    AIC = 2*p + negloglik
    #AIC = 2 * p + sum(ni) * log(parameters$sigma2)
    
    return(AIC)
}
