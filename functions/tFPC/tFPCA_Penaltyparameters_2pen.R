simplex <- function(K_fold,z,point,edges,n,d,r,J,p,station,ni,basismatrix,Pt,init,char){
  source("../functions/tFPC/CV_tFPCA.R")
  alpha=0.5 #default 0.5 0.5 0.5
  beta=0.5
  gamma=0.5
  
  design = matrix(ncol=2,nrow=3)
  design[1,]=c(init[1],0.0001) #default c(5,0) #Jun 3, 2020
  design[2,]=c(0.0001,init[2]) #default c(0,5) #Jun 3, 2020
  design[3,]=0.618*c(init[1],init[2])
  y=rep(0,3)
  for(i in 1:3){
    y[i] = CrossVal(K_fold,z,point,edges,n,d,r,J,p,station,ni,basismatrix,Pt,design[i,1],0,design[i,2],char)
  }
  
  
  ## Simplex Algorithm
  iter=1
  print("the criteria is")
  print(sd(y/(mean(y)+0.00001)))

  while(sd(y/(mean(y)+0.00001))>0.0001 & iter<20) #0.05 change to 0.01.
  {
    print(paste("start the ", iter, "-th iteration",sep=""))
    # ind_h=(1:3)[y==max(y)]
    # ind_l=(1:3)[y==min(y)]
    ord = order(y)
    ind_l=ord[1]
    ind_h=ord[3]
    barP = apply(design[-ind_h,],2,mean)
    barP[barP<0]=0
    starP = (1+alpha)*barP-alpha*design[ind_h,]
    starP[starP<0]=0
    
    print(paste("currently run", starP[1],starP[2], sep=" "))
    stary=CrossVal(K_fold,z,point,edges,n,d,r,J,p,station,ni,basismatrix,Pt,starP[1],0,starP[2],char)
    if(stary<y[ind_l]){
      ssP=(1+gamma)*starP - gamma*barP
      ssP[ssP<0]=0
      print(paste("currently run", ssP[1], ssP[2], sep=" "))
      ssy=CrossVal(K_fold,z,point,edges,n,d,r,J,p,station,ni,basismatrix,Pt,ssP[1],0,ssP[2],char)
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
        print(paste("currently run", ssP[1],ssP[2], sep=" "))
        ssy=CrossVal(K_fold,z,point,edges,n,d,r,J,p,station,ni,basismatrix,Pt,ssP[1],0,ssP[2],char)
        if(ssy>y[ind_h]){
          design=(design+design[ind_l,])/2
          for(i in 1:3)
          {
          	print(paste("currently run", design[i,1],design[i,2], sep=" "))
            y[i]=CrossVal(K_fold,z,point,edges,n,d,r,J,p,station,ni,basismatrix,Pt,design[i,1],0,design[i,2],char)
          }
        }else{
          design[ind_h,]=ssP
          y[ind_h]=ssy
        }
      }
    }
    #print(cbind(design,y))

    iter=iter+1
    print("the criteria is")
  	print(sd(y/(mean(y)+0.00001)))
  }
  index=which(y==min(y))
  penalty = design[index[1],]
  result = list(penalty, iter)
  return(result)
}
