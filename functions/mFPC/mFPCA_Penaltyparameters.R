simplex <- function(K_fold,J,z,rawbasis,ni,Q2,K){
  source("../functions/mFPC/CV_mFPCA.R")
  alpha=0.5  #default 1
  beta=0.5
  gamma=0.5 # default 2
  
  design = matrix(ncol=2,nrow=3)
  design[1,]=c(5,0.01) #default c(5,0) #Jun 3, 2020
  design[2,]=c(0.01,5) #default c(0,5) #Jun 3, 2020
  design[3,]=0.618*c(5,5)
  y=rep(0,3)
  for(i in 1:3){
    y[i] = CrossVal(K_fold,J,z,rawbasis,ni,Q2,K,design[i,1],design[i,2])
  }
  
  ## Simplex Algorithm
  iter=1
  while(sd(y/(mean(y)+0.0001))>0.05 & iter<20)
  {
    # ind_h=(1:3)[y==max(y)]
    # ind_l=(1:3)[y==min(y)]
    ord = order(y)
    ind_l=ord[1]
    ind_h=ord[3]
    barP = apply(design[-ind_h,],2,mean)
    barP[barP<0]=0
    starP = (1+alpha)*barP-alpha*design[ind_h,]
    starP[starP<0]=0
    stary=CrossVal(K_fold,J,z,rawbasis,ni,Q2,K,starP[1],starP[2])
    if(stary<y[ind_l]){
      ssP=(1+gamma)*starP
      ssP[ssP<0]=0
      ssy=CrossVal(K_fold,J,z,rawbasis,ni,Q2,K,ssP[1],ssP[2])
      if(ssy<y[ind_l]){
        design[ind_h,]=ssP
        y[ind_h]=ssy
      }else{
        design[ind_h,]=starP
        y[ind_h]=stary
      }
    }else{
      if(stary<=y[ord[2]]){
        design[ind_h,]=starP
        y[ind_h]=stary
      }else{
        if(stary<=y[ind_h]){
          design[ind_h,]=starP
          y[ind_h]=stary
        }
        ssP=beta*design[ind_h,]+(1-beta)*barP
        ssP[ssP<0]=0
        ssy=CrossVal(K_fold,J,z,rawbasis,ni,Q2,K,ssP[1],ssP[2])
        if(ssy>y[ind_h]){
          design=(design+design[ind_l,])/2
          for(i in 1:3)
          {
            y[i]=CrossVal(K_fold,J,z,rawbasis,ni,Q2,K,design[i,1],design[i,2])
          }
        }else{
          design[ind_h,]=ssP
          y[ind_h]=ssy
        }
      }
    }
    #print(cbind(design,y))
    iter=iter+1
  }
  index=which(y==min(y))
  penalty = design[index[1],]
  result = list(penalty, iter)
  return(result)
}
