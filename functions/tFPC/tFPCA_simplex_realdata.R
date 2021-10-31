simplex <- function(K_fold,basis,P,basismatrix,Pt,z,point,edges,n,d,r,J,p,ni,init,char){
  source("../functions/tFPC/CV_tFPCA_realdata.R")
  alpha=0.5
  beta=0.5
  gamma=0.5
  
  design = matrix(ncol=3,nrow=4)
  design[1,]=c(init[1],0.000001,init[3])
  design[2,]=c(init[1],init[2],0.0001)
  design[3,]=c(0.0001,init[2],init[3])
  design[4,]=0.618*c(init[1],init[2],init[3])
  y=rep(0,4)
  for(i in 1:4){
      print(c("initial",i))
      print(design[i,])
    y[i] = try(CrossVal(K_fold,basis,P,basismatrix,Pt,z,point,edges,n,d,r,J,p,ni,design[i,1],design[i,2],design[i,3],char))

    if(is.nan(y[i])){y[i] = 1e10}
  }
  
  print(cbind(design,y))
  
  print("the criteria is")
  print(sd(y/(mean(y)+0.00001)))
  
  ## Simplex Algorithm
  iter=1
  while(sd(y/(mean(y)+0.00001))>0.05 & iter<20) #change the convergence condition
  {
    print(paste("start the ", iter, "-th iteration",sep=""))
    
    # ind_h=(1:3)[y==max(y)]
    # ind_l=(1:3)[y==min(y)]
    ord = order(y)
    ind_l=ord[1]
    ind_h=ord[4]
    
    #### testing #####
    print(design[ind_l,])
    
    barP = apply(design[-ind_h,],2,mean)
    barP[barP<0]=0
    starP = (1+alpha)*barP-alpha*design[ind_h,]
    starP[starP<0]=0
    print("start stary")
    print(starP)
    stary=try(CrossVal(K_fold,basis,P,basismatrix,Pt,z,point,edges,n,d,r,J,p,ni,starP[1],starP[2],starP[3],char))
    
    if(is.nan(stary)){stary = 1e10}
    
    
    #test
    print("starP")
    print(starP)
    
    if(stary<y[ind_l]){
      print("stary < y[ind_l]")
      ssP=(1+gamma)*starP - gamma*barP
      ssP[ssP<0]=0
      print("start ssy")
      print("ssP")
      print(ssP)
      ssy=try(CrossVal(K_fold,basis,P,basismatrix,Pt,z,point,edges,n,d,r,J,p,ni,ssP[1],ssP[2],ssP[3],char))
      if(is.nan(ssy)){ssy = 1e10}
      
      #test
      
      if(ssy<y[ind_l]){
        print("ssy < y[ind_l]")
        design[ind_h,]=ssP
        y[ind_h]=ssy
      }else{
        design[ind_h,]=starP
        y[ind_h]=stary
      }
    }else{
      if(stary<=y[ord[2]]){
        print("stary <= y[ord[2]]")
        design[ind_h,]=starP
        y[ind_h]=stary
      }else{
        if(stary<=y[ind_h]){
          design[ind_h,]=starP
          y[ind_h]=stary
        }
        ssP=beta*design[ind_h,]+(1-beta)*barP
        ssP[ssP<0]=0
        ssy=try(CrossVal(K_fold,basis,P,basismatrix,Pt,z,point,edges,n,d,r,J,p,ni,ssP[1],ssP[2],ssP[3],char))
        if(is.nan(ssy)){ssy = 1e10}
        
        #test
        print("ssy")
        print(ssy)
        
        if(ssy>y[ind_h]){
          design=(design+design[ind_l,])/2
          for(i in 1:4)
          {
            y[i]=try(CrossVal(K_fold,basis,P,basismatrix,Pt,z,point,edges,n,d,r,J,p,ni,design[i,1],design[i,2],design[i,3],char))
            if(is.nan(y[i])){y[i] = 1e10}
          }
        }else{
          design[ind_h,]=ssP
          y[ind_h]=ssy
        }
      }
    }

    print(cbind(design,y))

    print("the criteria is")
    print(sd(y/(mean(y)+0.00001)))
    iter=iter+1
  }
  index=which(y==min(y))
  penalty = design[index[1],]
  print("final penalty")
  print(penalty)
  result = list(penalty, iter)
  return(result)
}


##cannot improve the computational speed anymore, this is the best.
