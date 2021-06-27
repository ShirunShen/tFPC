#### EM algorithm for mFPCA ####

## initialize the value ##
EMinitial <- function(z,ni,basis,J){
  lambda = 0.001 #default 1  ##ridege parameter
  f = {}
  sigma20 = 0
  for(i in 1:length(ni)){
    Bi = basis[indfunc(i,ni),]
    zi = z[indfunc(i,ni)]
    fi = solve(t(Bi)%*%Bi+lambda*diag(ncol(Bi)))%*%t(Bi)%*%zi
    sigma20 = sigma20 + sum((zi - Bi%*%fi)^2)/(ni[i])
    f  = cbind(f,fi)
    
  }
  
  fbar = apply(f,1,mean)  ##average of parameters ridge regression 
  fcen = f                ##centralization for f matrix
  for(i in 1:ncol(f)){
    fcen[,i] = f[,i] - fbar
  }
  
  
  ##initial value of sigma2,thetamu,Theta,D
  sigma20 = sigma20/length(ni)
  #thetamu0 = as.matrix(fbar)
  thetamu0 = solve(t(basis) %*% basis) %*% t(basis) %*% z #Jun 10, 2021
  sv = svd(fcen,nu=J)      #rank-J reduced SVD
  Theta0 = sv$u
  D0 = diag(sv$d[1:J]^2)/length(ni)
  
  return(list(sigma20 = sigma20, thetamu0 = thetamu0, Theta0 = Theta0, D0 = D0))
}
