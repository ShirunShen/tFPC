simplex <- function(K_fold,J,z,basis,P,ni,init){
  source("../functions/mFPC/CV_mFPCA_realdata.R")
  alpha=1
  beta=0.5
  gamma=2
  
  design = matrix(ncol=2,nrow=3)
  design[1,]=c(init[1],0.1)
  design[2,]=c(0.1,init[2])
  design[3,]=0.618*c(init[1],init[2])
  y=rep(0,3)
  for(i in 1:3){
    y[i] = CrossVal(K_fold,J,z,basis,P,ni,design[i,1],design[i,2])
  }
  
  
  
  ## Simplex Algorithm
  iter=1
  while(sd(y/(mean(y)+0.01))>0.05 & iter<50)
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
    stary=CrossVal(K_fold,J,z,basis,P,ni,starP[1],starP[2])
    if(stary<y[ind_l]){
      ssP=(1+gamma)*starP
      ssP[ssP<0]=0
      ssy=CrossVal(K_fold,J,z,basis,P,ni,ssP[1],ssP[2])
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
        ssy=CrossVal(K_fold,J,z,basis,P,ni,ssP[1],ssP[2])
        if(ssy>y[ind_h]){
          design=(design+design[ind_l,])/2
          for(i in 1:3)
          {
            y[i]=CrossVal(K_fold,J,z,basis,P,ni,design[i,1],design[i,2])
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
