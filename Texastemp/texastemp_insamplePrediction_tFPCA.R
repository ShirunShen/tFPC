rm(list=ls())
setwd("~/Desktop/tFPC_final/Texastemp/")

library(Matrix)
library(MASS)
library(mgcv)
library(ape)
library(fields)
library(fda)

source("../functions/tFPC/Multi_tFPCA.R")

##### loading data
temperature = as.matrix(read.table(file="./Data/temperature.txt")) #the last column is the yearly average temperature
stationid = read.table(file="./Data/stationid.txt")$x
years = read.table(file="./Data/years.txt")$x
stationlocations = read.table(file="./Data/stationlocation.txt")
x = stationlocations[,3]
y = stationlocations[,2]


#usaboundary = read.table(file="./Data/usaboundary.txt")
#pdf("usa.pdf")
#plot(usaboundary,type="l")
#points(x,y,col="red")
#dev.off()
#stationlocations[1013,]


# if change the region of Texas, from here.
#Texasboundary=list(read.table("./Data/texasboundary.dat",header=T))
Texasboundary=list(read.table("./Data/texasregion.txt",header=T))
names(Texasboundary[[1]])=c("x","y")

ind_texas = which(inSide(Texasboundary,x,y))



#ind_texas = which(stationlocations[,4]=="TX")
Texasstation = stationlocations[ind_texas,3:2]

id_texas = stationlocations[ind_texas,1]
##### pre-processing data to make it as what we expected.

temp_texas = NULL
year_texas = NULL
num_year   = NULL
for(i in 1:length(id_texas)){
    ind = grep(id_texas[i],stationid)
    ind_year = years[ind]
    year_texas = c(year_texas,ind_year)
    temp_texas = rbind(temp_texas,temperature[ind,1:12])
    num_year = c(num_year,length(ind))
}
ind_station = NULL
for(i in 1:length(num_year)){
    ind_station = c(ind_station,rep(i,num_year[i]))
}
temp_texas = cbind(temp_texas,year_texas,ind_station)



texastemp1 = matrix(-9999,ncol=length(ind_texas),nrow=(2014-1878+1)*12)
for(i in 1:nrow(temp_texas)){
    ind_station = temp_texas[i,14]
    year = temp_texas[i,13]
    #print(year)
    texastemp1[(year-1878)*12 + 1:12, ind_station] = as.vector(temp_texas[i,1:12])
}

for(i in 1:length(texastemp1)){
    if(texastemp1[i] < -100){
        texastemp1[i] = NA
    }
}

CelsiusTrans <- function(x){
    return((x-32)/1.8)
}
texastemp1 = CelsiusTrans(texastemp1/10)

Texastemp = texastemp1[(37*12+1):nrow(texastemp1),] #start from 1915 to 2014, total 100years


#pdf("test.pdf")
#for(i in 1:ncol(Texastemp)){
#    plot(Texastemp[1:38,i],type="o",xaxt="n")
#    axis(1, at=c(13,25,37), labels=c(13,25,37))
#}
#dev.off()


#pdf("heatmap_data_texastemp.pdf")
#fields::image.plot(texastemp)
#dev.off()

##### plot triangulation of Texas maps

V=matrix(c(-106.74719, 31.99396,
-103.07511, 31.99396,
-103.04646, 36.49,
-100.00405,36.5,
-99.99833, 34.56081,
-96.91581, 33.96494,
-94.0967, 33.98678,
-93.30778, 30.05110,
-97.23559, 27.79319,
-97.09073, 25.7812,
-99.63558, 27.45041,
-101.63385, 29.78235,
-103.66106, 28.58594,
-99.80000, 31.5000,
-97.00000, 30.05000,
-99.20000, 26.30000),ncol=2,byrow=T)

Tri=matrix(c(2,1,13,3,2,5,4,3,5,5,2,14,5,14,6,6,14,15,6,15,7,7,15,8,8,15,9,9,15,11,10,9,11,16,10,11,11,15,12,12,15,14,12,14,2,12,2,13),ncol=3,byrow=T)



# 38 NEW BRUNFELS

Texasstation38 = Texasstation[38,]
Texasstation_38 = Texasstation[-38,]

Texastemp38 = Texastemp[,38]
Texastemp_38 = Texastemp[,-38]
Texastemp = Texastemp_38

ind_station = NULL

#z = NULL
#station = NULL
#ni = NULL
#for(i in 1:nrow(Texastemp)){
#    count = 0
#    for(j in 1:ncol(Texastemp)){
#        if(!is.na(Texastemp[i,j])){
#            z = c(z,Texastemp[i,j])
#            count = count + 1
#            station = rbind(station,Texasstation_38[j,])

#            ind_station = c(ind_station,j)
#        }
#    }
#    ni = c(ni,count)
#} # need some time!!!

#write.table(z, file="./Data/temp_38.txt")
#write.table(ni, file="./Data/ni_38.txt")
#write.table(station, file="./Data/station_38.txt")

z = as.matrix(read.table("./Data/temp_38.txt"))
ni = as.matrix(read.table("./Data/ni_38.txt"))
station = as.matrix(read.table("./Data/station_38.txt"))


z_tFPCA = z

######## area of Texas #####
m=70;n=70
xm=seq(min(V[,1]),max(V[,1]),length=m)
yn=seq(min(V[,2]),max(V[,2]),length=n)
x=rep(xm,n)
y=rep(yn,rep(m,n))
ind=inSide(Texasboundary,x,y) ##index of (xx,yy) inside Texas, understood!

#calculate area of texas
#area=emp.area(Texasboundary,10^6)

area = 64.91116

d=3
r=1

H=smoothness(V,Tri,d,r)
qrdecom=qr(t(H))
R=qr.R(qrdecom)
Q=qr.Q(qrdecom,TRUE)
nh=nrow(H)
Q2=Q[,(qrdecom$rank+1):ncol(Q)] #Q2 constraint


B1 = beval(V,Tri,d,x[ind],y[ind])
R = qr.R(qr(B1%*%Q2))
omat = solve(R)*sqrt(length(x[ind])/area)
omat = Q2%*%omat #Gram-Schmidt Orthonormalization



Bstation = beval(V,Tri,d,Texasstation[,1],Texasstation[,2])
rawbasis2 = Bstation%*%omat

rawbasis238 = rawbasis2[38,]
rawbasis2_38 = rawbasis2[-38,]
rawbasis2 = rawbasis2_38

#basis = NULL
#for(i in 1:nrow(Texastemp)){
#    for(j in 1:ncol(Texastemp)){
#        if(!is.na(Texastemp[i,j])){
#            basis = rbind(basis,rawbasis2_38[j,])
#        }
#    }
#}

#write.table(basis,file="./Data/basis_38.txt")

basis = as.matrix(read.table("./Data/basis_38.txt"))




#save(file="inSample_prediction.RData")

### continue ###
#setwd("~/Desktop/TemporalFPCAproject/NewTexastemp/")

#library(Matrix)
#library(MASS)
#library(mgcv)
#library(ape)
#library(fields)
#library(fda)

#load(file="inSample_prediction.RData")



#penalty matrix
Ene=energy(V,Tri,d,1)
P=t(omat)%*%Ene%*%(omat)  #with orthogonalizations



####### tFPCA #######

set.seed(1)
J=3    # number of principal component vectors

splinebasis = create.fourier.basis(c(0,length(ni)), nbasis=51, period = 12)
basismatrix = eval.basis(0:(length(ni)-1), splinebasis)


#test, April 5
basismatrix = cbind(basismatrix,1:nrow(basismatrix),(1:nrow(basismatrix))^2,(1:nrow(basismatrix))^3)

test = t(basismatrix) %*% basismatrix


names(basismatrix) = NULL
Pt = fourierpen(splinebasis)

#test, April 5
#Pt = cbind(Pt,0)
#Pt = rbind(Pt,0)
Pt = cbind(Pt,0,0,0)
Pt = rbind(Pt,0,0,0)


# qr decomposition to make it orthonormal
#basisprime = qr(basismatrix)
#basismatrix1 = qr.Q(basisprime)
#orthoterm = qr.R(basisprime)

#Pt1 = t(solve(orthoterm)) %*% Pt %*% solve(orthoterm)

#basismatrix = basismatrix1
#Pt = Pt1




nb = ncol(basis)
nc = ncol(basismatrix)
I_nb = diag(nb)
p=2


#### select p ####
source("../functions/tFPC/AIC.R")

set.seed(2)
#AIC_tFPCA = NULL
#for(p in 1:3){
#    temp = AIC(p,K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,ni)
#    AIC_tFPCA = c(AIC_tFPCA, temp)
#}
#p = which(AIC_tFPCA == min(AIC_tFPCA))


### test, Oct 26, 2021

source("../functions/tFPC/CV_tFPCA_realdata.R")
grid_mus = c(1,100,10000)
grid_mut = c(1,100,10000)
grid_pc = c(1,100,10000)

grid_mupc = expand.grid(grid_mus, grid_mut, grid_pc)


one_CV <- function(i,basis,P,basismatrix,Pt,z,V,Tri,t,d,r,J,p,ni,grid_mupc){
    tmp=try({CrossVal(K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,t,d,r,J,p,ni,grid_mupc[i,1],grid_mupc[i,2],grid_mupc[i,3],"EMalgorithm")},silent=TRUE)
    tmp = ifelse(class(tmp)=="try-error",10000,tmp)
    print(i)
    print(c(tmp,unlist(grid_mupc[i,])))
    return(c(tmp,unlist(grid_mupc[i,])))
}


penalty_result = NULL
for(i in 1:nrow(grid_mupc)){
    tmp = one_CV(i,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,p,ni,grid_mupc)
    penalty_result <- rbind(penalty_result, tmp)
}

write.table(penalty_result,file="insamplePrediction_penalty.txt")
init = penalty_result[which.min(penalty_result[,1]),-1]

#init = c(10000,10000,100)

### end test, Oct 26, 2021



set.seed(2)
#init = c(1,10000,10000) #(1,1,5000)
source("../functions/tFPC/tFPCA_simplex_realdata.R")
penalty = simplex(K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,p,ni,init,"EMalgorithm")
#save.image("2.RData")
#penalty = list(c(6180,0.0000618,6180)) #June 8, 2021
#penalty = list(c(0.618,0.618,1236),1) #default 0.1, 1, 10000
#penalty = list(c(0.618/(1526^2), 0.618*(1526^2), 1236),1)

#penalty = list(c(0.1545,4847.9,6545.077),1) 
#penalty = list(c(0.809,3090,8090),1)
#penalty = list(c(1,0.00008090,13090),1) #Best one.
lambmus = penalty[[1]][1]
lambmut = penalty[[1]][2]
lambpc  = penalty[[1]][3]


#source("../functions/tFPC/CV_tFPCA_realdata.R")
#set.seed(2)
#CrossVal(K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,p,ni,lambmus,lambmut,lambpc,"EMalgorithm")


set.seed(5)   #default 6, # 3, 6
K0 = rep(0.5,p)
sigma20 = 1.0
HJ0 = diag(J)
Theta0 = matrix(runif(J*nb,-1,1),ncol=J)
thetab0 = runif(nb,-1,1)
thetab0 = thetab0/norm(thetab0)
thetac0 = runif(nc,-1,1)
thetac0 = thetac0/norm(thetac0)

Rcpp::sourceCpp("../functions/tFPC/EMalgorithm.cpp")

#### for testing the initial values ####
### based on EMalgorithm_tFPCA0_manifoldthetab.cpp
set.seed(1)
lambdab = 0.0001 #default = 1
lambdac = 0.0001
test = EMinit(z,basis,basismatrix,ni,thetab0,thetac0,P,Pt,lambdab,lambdac)
thetab0 = test$thetab
thetac0 = test$thetac


mu0 = basis[1:ni[1],] %*% thetab0 %*% t(thetac0) %*% basismatrix[1,]
for(t in 2:length(ni)){
    temp = basis[sum(ni[1:(t-1)]) + 1:ni[t], ] %*% thetab0 %*% t(thetac0) %*% basismatrix[t,]
    mu0 = c(mu0,temp)
}
pcpart = z - mu0

f = NULL
for(t in 1:length(ni)){
    Bt = basis[indfunc(t,ni),]
    pcpart_t = pcpart[indfunc(t,ni)]
    ft = solve(t(Bt) %*% Bt +  0.001 * diag(ncol(Bt))) %*% t(Bt) %*% pcpart_t
    f = cbind(f,ft)
}

sv = svd(f,nu=J)
Theta0 = sv$u
HJ0 = diag(sv$d[1:J]^2)/length(ni)


mutest = basis[1:ni[1],] %*% thetab0 %*% t(thetac0) %*% basismatrix[1,]
for(t in 2:length(ni)){
    temp = basis[sum(ni[1:(t-1)]) + 1:ni[t], ] %*% thetab0 %*% t(thetac0) %*% basismatrix[t,]
    mutest = c(mutest, temp)
}
residual = z - mutest;
in_sample_tMSE = sum(residual^2)/length(residual)
in_sample_tMAE = sum(abs(residual))/length(residual)
in_sample_tMAE #manifold thetac = 6.82, manifold thetab = 16.15
in_sample_tMSE #manifold thetac = 68.21, manifold thetab = 320.45




system.time({parameters = EMalgorithm(z,basis,basismatrix,P,Pt,ni,HJ0,thetab0,thetac0,
    Theta0,K0,sigma20,lambmus,lambmut,lambpc)})

thetabhat = parameters$thetab
thetachat = parameters$thetac
Thetahat = parameters$Theta

####### muhat is ready to visualize #########
muhat = basis[1:ni[1],] %*% thetabhat %*% t(thetachat) %*% basismatrix[1,]

for(t in 2:length(ni)){
    temp = basis[sum(ni[1:(t-1)]) + 1:ni[t], ] %*% thetabhat %*% t(thetachat) %*% basismatrix[t,]
    muhat = c(muhat, temp)
    
}

######## zhat is ready to visualize #########
temp = rep(0,length(muhat))
for(j in 1:J){
    alphajhat = parameters$alpha[seq(j,J*length(ni),by=J)]
    temp = temp + (basis %*% parameters$Theta[,j]) * rep(alphajhat,ni)
}

zhat = muhat + temp

####### tMISE #######
residual = z - zhat;
in_sample_tMSE = sum(residual^2)/length(residual)
in_sample_tMAE = sum(abs(residual))/length(residual)

####### The whole plot ########

basis1 = B1 %*% omat  #B1 = beval(V,Tri,d,x[ind],y[ind]) # for every time point

#test = matrix(basis1, nrow=nrow(basis1)*length(ni),ncol=ncol(basis))

muhat1 = NULL
for(t in 1:length(ni)){
    muhat1 = c(muhat1, basis1 %*% thetabhat %*% t(thetachat) %*% basismatrix[t,])
}

temp = rep(0,length(muhat1))
for(j in 1:J){
    alphajhat = parameters$alpha[seq(j,J*length(ni),by=J)]
    temp = temp +  rep(basis1 %*% parameters$Theta[,j],length(ni)) * rep(alphajhat,rep(nrow(basis1),length(ni)))
}
zhat1 = muhat1 + temp

nnp = nrow(basis1)

ym = rep(NA,m*n)

zhat1_tFPCA = zhat1
summary(zhat1_tFPCA)


muhat38 = rawbasis238 %*% thetabhat * t(thetachat) %*% basismatrix[1,]
for(t in 2:length(ni)){
    muhat38 = c(muhat38, rawbasis238 %*% thetabhat * t(thetachat) %*% basismatrix[t,])
}

temp = rep(0, length(muhat38))
for(j in 1:J){
    alphajhat = parameters$alpha[seq(j,J*length(ni),by=J)]
    temp = temp + as.numeric(rawbasis238 %*% Thetahat[,j]) * alphajhat
}
zhat38 = muhat38 + temp

zhat38_tFPCA = zhat38

Pred_error_tFPCA = (zhat38_tFPCA - as.vector(Texastemp38))
PPE_tFPCA = (zhat38_tFPCA - as.vector(Texastemp38))/as.vector(Texastemp38)
summary(Pred_error_tFPCA)


summary(zhat38)

#save.image("2.RData")


MSE_tFPCA = NULL
MAE_tFPCA = NULL

for(ind in 1:12){
    singlemonth = seq(ind,nrow(Texastemp),by=12)
    MSE_tFPCA = c(MSE_tFPCA, mean(Pred_error_tFPCA[singlemonth]^2,na.rm=TRUE))
    MAE_tFPCA = c(MAE_tFPCA, mean(abs(Pred_error_tFPCA[singlemonth]),na.rm=TRUE))
}
MAE_tFPCA



tMSE_tFPCA = mean((Texastemp38-zhat38_tFPCA)^2,na.rm=TRUE)
tMSE_tFPCA



setwd("~/Desktop/tFPC_final/Texastemp/")

library(Matrix)
library(MASS)
library(mgcv)
library(ape)
library(fields)
library(fda)
load(file="inSample_prediction.RData")


MSE_tFPCA = NULL
MAE_tFPCA = NULL

for(ind in 1:12){
    singlemonth = seq(ind,nrow(Texastemp),by=12)
    MSE_tFPCA = c(MSE_tFPCA, mean(Pred_error_tFPCA[singlemonth]^2,na.rm=TRUE))
    MAE_tFPCA = c(MAE_tFPCA, mean(abs(Pred_error_tFPCA[singlemonth]),na.rm=TRUE))
}


MSE_mFPCA = read.table(file="z38mse_mFPCA.txt")
MAE_mFPCA = read.table(file="z38mae_mFPCA.txt")

MSE_tFPCA - MSE_mFPCA
MAE_tFPCA - MAE_mFPCA


#write.table(rbind(MeanSE_tFPCA,MeanSE_mFPCA),"Prediction_realdata.txt")

pdf("Prediction_realdata.pdf",width=8,height=4)
plot(MAE_tFPCA, type="b",ylim=c(0.4,1.0),pch=19,xlab="Month",ylab="MAE",col="red",xaxt="n")
axis(1, at = 1:12,label=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
lines(MAE_mFPCA, pch=18, type="b", col="blue",lty=2)
legend(10,1,legend=c("tFPC","mFPC"),col=c("red","blue"),lty=1:2, cex=0.9)
dev.off()





save.image(file="inSample_prediction.RData")


# 
# 
# set.seed(1)
# zhat38_mFPCA_summary = NULL
# PPE_mFPCA_summary = NULL
# MSE_mFPCA = NULL
# MAE_mFPCA = NULL
# 
# 
# for(ind in 1:12){
#     print(ind)
#     singlemonth = seq(ind,nrow(Texastemp),by=12)
#     Texastemp1 = Texastemp[singlemonth,] #single month
#     
#     
#     ind_singlemonth = NULL
#     
#     if(singlemonth[1] == 1){
#         ind_singlemonth = 1:ni[1]
#     }else{
#         ind_singlemonth = sum(ni[1:(singlemonth[1]-1)])+1:ni[singlemonth[1]]
#     }
#     
#     
#     #pdf("rawdata.pdf")
#     #image(as.matrix(Texastemp1))
#     #dev.off()
#     
#     
#     ind_station = NULL
#     
#     z = NULL
#     station = NULL
#     ni = NULL
#     for(i in 1:nrow(Texastemp1)){
#         count = 0
#         for(j in 1:ncol(Texastemp1)){
#             if(!is.na(Texastemp1[i,j])){
#                 z = c(z,Texastemp1[i,j])
#                 count = count + 1
#                 station = rbind(station,Texasstation_38[j,])
#                 
#                 ind_station = c(ind_station,j)
#             }
#         }
#         ni = c(ni,count)
#     } # need some time!!!
#     
#     z_mFPCA = z
#     
#     
#     basis = NULL
#     for(i in 1:nrow(Texastemp1)){
#         for(j in 1:ncol(Texastemp1)){
#             if(!is.na(Texastemp1[i,j])){
#                 basis = rbind(basis,rawbasis2[j,])
#             }
#         }
#     }
#     #if(basis[1,1] < 0){basis = -basis}
#     
#     
#     
#     
#     source("../functions/mFPC/EMinitial_mFPCA_realdata.R")
#     Rcpp::sourceCpp("../functions/mFPC/EMalg.cpp")
#     
#     source("../functions/mFPC/CV_mFPCA_realdata.R")
#     source("../functions/mFPC/mFPCA_simplex_realdata.R")
#     
#     init_mFPCA = c(10,10) #default 10 10
#     penaltyPar = simplex(K_fold=5,J,z,basis,P,ni,init_mFPCA)
#     
#     #source("../functions/modelbasedPCA/CV_mFPCA_realdata.R")
#     #CrossVal(K_fold=5,J,z,basis,P,ni,1,1)
#     
#     #penaltyPar = list(c(2,1))
#     
#     lambmu = penaltyPar[[1]][1]
#     lambpc = penaltyPar[[1]][2]
#     
#     ####
#     EMinit = EMinitial(z,ni,basis,J)
#     sigma20 = EMinit$sigma20
#     D0 = EMinit$D0
#     thetamu0 = EMinit$thetamu0
#     Theta0 = EMinit$Theta0
#     
#     EM = EM_algorithm(z,basis,P,ni,lambmu,lambpc,sigma20,D0,thetamu0,Theta0)
#     
#     thetamu_mFPCA = EM$thetamu
#     Thetahat_mFPCA = EM$BigTheta
#     
#     muhat38_mFPCA = rawbasis238 %*% thetamu_mFPCA
#     muhat38_mFPCA = rep(muhat38_mFPCA,length(ni))
#     
#     temp = rep(0,length(muhat38_mFPCA))
#     for(j in 1:J){
#         temp = temp + EM$alpha[j,] * as.numeric(rawbasis238 %*% Thetahat_mFPCA[,j])
#     }
#     zhat38_mFPCA = muhat38_mFPCA + temp
#     zhat38_mFPCA_summary = rbind(zhat38_mFPCA_summary,zhat38_mFPCA)
#     Pred_error_mFPCA = (zhat38_mFPCA - as.vector(Texastemp38[singlemonth]))
#     PPE_mFPCA = (zhat38_mFPCA - as.vector(Texastemp38[singlemonth]))/as.vector(Texastemp38[singlemonth])
#     PPE_mFPCA_summary = rbind(PPE_mFPCA_summary,PPE_mFPCA)
#     MSE_mFPCA = c(MSE_mFPCA, mean(Pred_error_mFPCA^2,na.rm=TRUE))
#     
#     MAE_mFPCA = c(MAE_mFPCA, mean(abs(Pred_error_mFPCA),na.rm=TRUE))
# }
# 
# write.table(PPE_mFPCA_summary,file="test_PPE_mFPCA_summary.txt")
# write.table(MAE_mFPCA, file="test_MAE_mFPCA.txt")
# PPE_mFPCA_summary = unlist(read.table(file="test_PPE_mFPCA_summary.txt"))
# 
# 
# 
# month = rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),200)
# month = factor(month, levels = month.abb)
# method = c(rep("tFPCA",1200),rep("mFPCA",1200))
# PPE = c(PPE_tFPCA,as.vector(PPE_mFPCA_summary)) #prediction proportion error
# data_PPE = data.frame(month,method,PPE)
# 
# pdf("Prediction_boxplot1.pdf",width=8,height=4)
# library(ggplot2)
# ggplot(data_PPE,aes(x=month,y=PPE,fill=method)) + geom_boxplot()
# dev.off()
# 







# pdf("inSample_residual_Jun3.pdf",height=5,width=10)
# par(mar=c(5,5,2,2))
# plot((zhat38_tFPCA - Texastemp38)/Texastemp38,ylim=c(-6,6),type="p",col="red",pch=2,cex=0.5,xlab="Month",ylab="Prediction Error",xaxt="n",cex.axis=1.2,cex.lab=1.5)
# abline(h=0,col="black")
# char = c("Jan.1915","Jan.1935","Jan.1955","Jan.1975","Jan.1995","Dec.2014")
# axis(1,c(seq(1,1200,by=240),1200),char,cex.axis=1.2)#seq(0,1200,by=120))
# points((zhat38_mFPCA - Texastemp38)/Texastemp38,pch=5,col="blue",cex=0.5)
# legend(0,6,legend=c("tFPCA","mFPCA"), col=c("red","blue"),pch=c(2,5),cex=0.8)
# dev.off()

#ind_empty = which(is.na(Texastemp38))
#MPAE_mFPCA = sum(abs(zhat38_mFPCA[-ind_empty]-Texastemp38[-ind_empty])/Texastemp38[-ind_empty])/length(Texastemp38[-ind_empty])
#MPAE_tFPCA = sum(abs(zhat38_tFPCA[-ind_empty]-Texastemp38[-ind_empty])/Texastemp38[-ind_empty])/length(Texastemp38[-ind_empty])
#c(MPAE_mFPCA,MPAE_tFPCA)

#test = cbind((zhat38_mFPCA-Texastemp38)/Texastemp38, (zhat38_tFPCA-Texastemp38)/Texastemp38)
#colnames(test) <- c("mFPCA", "tFPCA")
#pdf("boxplot_residual.pdf")
#boxplot(test)
#abline(h=0,col="black")
#dev.off()

