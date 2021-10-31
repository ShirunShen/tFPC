rm(list=ls())
setwd("~/Desktop/tFPC_final/Texastemp")

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

#stationlocations[1013,]

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


#### forecasting ####
n_forecast = 12# number of forecasting months




Texastemp = texastemp1[(37*12+1):(nrow(texastemp1)-n_forecast),] #start from 1915 to 2013, total 80 years

Texastemp_out = texastemp1[nrow(texastemp1)-(n_forecast-1):0,]

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


station = Texasstation
boundary = Texasboundary
#pdf("data_texas.pdf",height=4,width=4)
#par(mar=c(4,4,2,2))
#plot(boundary[[1]],type="l",xlab="",ylab="",xlim=c(-107,-93),cex.lab=1.3)
#mtext("Longitude", side=1, line=2)
#mtext("Latitude", side=2, line=2)
#points(station[,1],station[,2],pch=2,cex=.5)
#plot.t(V,Tri,1,"l")
#dev.off()


####### Further processing ######

#ind_station = NULL

#z = NULL
#station = NULL
#ni = NULL
#for(i in 1:nrow(Texastemp)){
#    count = 0
#    for(j in 1:ncol(Texastemp)){
#        if(!is.na(Texastemp[i,j])){
#            z = c(z,Texastemp[i,j])
#            count = count + 1
#            station = rbind(station,Texasstation[j,])
#
#            ind_station = c(ind_station,j)
#        }
#    }
#    ni = c(ni,count)
#} # need some time!!!

#z_tFPCA = z

#write.table(z,file="./Data/temp.txt")
#write.table(ni,file="./Data/ni.txt")
#write.table(station,file="./Data/station.txt")

z = as.matrix(read.table("./Data/temp49_100years.txt"))

ni = as.matrix(read.table("./Data/ni49_100years.txt"))
station = as.matrix(read.table("./Data/station49_100years.txt"))

z_tFPCA = z



### depart the in and out samples
ni_forecast = sum(ni[(length(ni)-n_forecast+1):length(ni)])
z_in = z[1:(length(z)-ni_forecast)]
z_out = z[length(z)-(ni_forecast-1):0]
z = z_in
ni_in = ni[1:(length(ni)-n_forecast),1]
ni_out = ni[length(ni)-(n_forecast-1):0,]
ni = ni_in
station_in = station[1:(nrow(station)-ni_forecast),]
station_out = station[nrow(station)-(ni_forecast-1):0,]
station = station_in

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
#rawbasis2 = B1%*%omat



#basis = NULL
#for(i in 1:nrow(Texastemp)){
#    for(j in 1:ncol(Texastemp)){
#        if(!is.na(Texastemp[i,j])){
#            basis = rbind(basis,rawbasis2[j,])
#        }
#    }
#}
#if(basis[1,1] < 0){basis = -basis}
#write.table(basis,file="../Data/basis.txt")

basis = as.matrix(read.table("./Data/basis49_100years.txt"))

basis_in = basis[1:(nrow(basis)-ni_forecast),]
basis_out = basis[nrow(basis)-(ni_forecast-1):0,]
basis = basis_in


#z = z[-(1:sum(ni[1:480]))]

#station = station[-(1:sum(ni[1:480])),]

#basis = basis[-(1:sum(ni[1:480])),]
#ni = ni[-(1:480)]



#penalty matrix
Ene=energy(V,Tri,d,1)
P=t(omat)%*%Ene%*%(omat)  #with orthogonalizations




set.seed(1)
J=3    # number of principal component vectors

splinebasis = create.fourier.basis(c(0,length(ni)), nbasis=51, period = 12)
basismatrix = eval.basis(0:(length(ni)-1), splinebasis)

#test, April 5 / June 22
#basismatrix = cbind(basismatrix,1:nrow(basismatrix),(1:nrow(basismatrix))^2,(1:nrow(basismatrix))^3)

polybasis = bs(1:length(ni),#knots=c(66)*12,
               degree=3, intercept=FALSE)

basismatrix = cbind(basismatrix,polybasis)


names(basismatrix) = NULL
Pt = fourierpen(splinebasis)

#test, April 5
#Pt = cbind(Pt,0,0,0)
#Pt = rbind(Pt,0,0,0)

Pt = cbind(Pt, matrix(0,ncol=ncol(basismatrix)-51,nrow=nrow(Pt)))
Pt = rbind(Pt, matrix(0,nrow=ncol(basismatrix)-51,ncol=ncol(Pt)))




# svd decomposition to make it orthonormal
#basisprime = svd(basismatrix)
#basismatrix1 = basisprime$u
#orthoterm = diag(basisprime$d) %*% t(basisprime$v)

#Pt1 = orthoterm %*% Pt %*% t(orthoterm)

#basismatrix = basismatrix1
#Pt = Pt1




nb = ncol(basis)
nc = ncol(basismatrix)
I_nb = diag(nb)


source("../functions/tFPC/AIC.R")

set.seed(2)
#AIC_tFPCA = NULL
#for(p in 1:3){
#    temp = AIC(p,K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,ni)
#    AIC_tFPCA = c(AIC_tFPCA, temp)
#}
#p = which(AIC_tFPCA == min(AIC_tFPCA))

p=2


### test Oct 2021
set.seed(1)
# source("../functions/tFPC/CV_tFPCA_realdata.R")
# grid_mus = c(1,100,10000)
# grid_mut = c(0.0001,0.01, 0.1)
# grid_pc = c(1,100,10000)
# 
# grid_mupc = expand.grid(grid_mus, grid_mut, grid_pc)
# 
# 
# one_CV <- function(i,basis,P,basismatrix,Pt,z,V,Tri,t,d,r,J,p,ni,grid_mupc){
#     tmp=try({CrossVal(K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,t,d,r,J,p,ni,grid_mupc[i,1],grid_mupc[i,2],grid_mupc[i,3],"EMalgorithm")},silent=TRUE)
#     tmp = ifelse(class(tmp)=="try-error",10000,tmp)
#     print(i)
#     print(c(tmp,unlist(grid_mupc[i,])))
#     return(c(tmp,unlist(grid_mupc[i,])))
# }
# 
# 
# penalty_result = NULL
# for(i in 1:nrow(grid_mupc)){
#     tmp = one_CV(i,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,p,ni,grid_mupc)
#     penalty_result <- rbind(penalty_result, tmp)
# }
# 
# init = penalty_result[which.min(penalty_result[,1]),-1]

#write.table(penalty_result,file="forecast_tFPCA_penalty.txt")

init = c(100,0.1,1)


### end test Oct 2021







set.seed(2)
#init = c(5000,0.0001,5000) #(1,1,5000)

source("../functions/tFPC/tFPCA_simplex_realdata.R")
penalty = simplex(K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,p,ni,init,"EMalgorithm")
#save.image("3.RData")
#penalty = list(c(6960,0.0001191,6960),3)


lambmus = penalty[[1]][1]
lambmut = penalty[[1]][2]
lambpc  = penalty[[1]][3]


#source("../functions/TemporalPCA/CV_tFPCA_realdata.R")
#set.seed(2)
#CrossVal(K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,p,ni,lambmus,lambmut,lambpc,"EMalgorithm_tFPCA0_manifoldthetab")



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
lambdab = 0.0001 #default = 10000
lambdac = 0.0001
test = EMinit(z,basis,basismatrix,ni,thetab0,thetac0,P,Pt,lambdab,lambdac)
thetab0 = test$thetab
thetac0 = test$thetac





system.time({parameters = EMalgorithm(z,basis,basismatrix,P,Pt,ni,HJ0,thetab0,thetac0,
    Theta0,K0,sigma20,lambmus,lambmut,lambpc)})

thetabhat = parameters$thetab
thetachat = parameters$thetac
Thetahat = parameters$Theta

#plot(basismatrix %*% thetachat,type="l")
#dev.off()


#### forecasting, April 9, 2020 ######
basis_T12 = eval.basis(length(ni)+0:(n_forecast-1),splinebasis)

polybasis_prime = predict(polybasis, newx=length(ni)+1:12)

basis_T12 = cbind(basis_T12, polybasis_prime)




# with trend, April, 2020, 1:948/960,(1:948)^2/(960^2))
#basis_T12 = cbind(basis_T12, (nrow(basismatrix)+1):(nrow(basismatrix)+n_forecast), ((nrow(basismatrix)+1):(nrow(basismatrix)+n_forecast))^2, ((nrow(basismatrix)+1):(nrow(basismatrix) + n_forecast))^3)


mu2hat_T12_pred = basis_T12 %*% thetachat

alpha_in12 = parameters$alpha

for(t in 1:n_forecast){
    temp = rep(0,J)
    for(i in 1:p){
        alpha_temp = alpha_in12[length(alpha_in12)-J*(i)+1:J]
        #print(alpha_temp)
        temp = temp + (parameters$K[i] * alpha_temp)
        #print(temp)
    }
    alpha_in12 = c(alpha_in12,temp)
}

alpha_T12 = matrix(alpha_in12[-(1:(length(alpha_in12)-J*n_forecast))],ncol=n_forecast)

muhat_pred = rawbasis2 %*% thetabhat %*% t(mu2hat_T12_pred)

z_pred = muhat_pred ### temporary
for(t in 1:n_forecast){
    for(j in 1:J){
        z_pred[,t] = z_pred[,t] + (rawbasis2 %*% parameters$Theta[,j]) * alpha_T12[j,t]
    }
}

temp_forecast12 = as.matrix(z_pred)

temp_true12 = t(Texastemp_out)

write.table(temp_forecast12, file="tFPCA2014.txt")

save.image(file="forecast.RData")

#mean_monthly = apply(temp_forecast12 - temp_true12, 2, function(x){return(mean(x,na.rm=TRUE))})
# J 0.74, F 1.20, M 1.26, A 0.30, M 0.89, J -0.41, J 0.25, A -0.44
# S 0.10, O -1.08, N 3.305, D -0.20





### plot ####

setwd("~/Desktop/tFPC_final/Texastemp")

temp_true = read.table(file="true2014.txt")
temp_mFPCA = read.table(file="mFPCA2014.txt")
temp_tFPCA = read.table(file="tFPCA2014.txt")

diff_mFPCA = unlist(temp_mFPCA) - unlist(temp_true)
diff_tFPCA = unlist(temp_tFPCA) - unlist(temp_true)

month = rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),each = 49)
month = c(month,month)
month = factor(month, levels = month.abb)
method = c(rep("mFPCA",49*12),rep("tFPCA",49*12))
PE = c(diff_mFPCA,diff_tFPCA) #prediction proportion error
data_diff = data.frame(month,method,PE)

pdf("Forecast_boxplot.pdf",width=8,height=4)
library(ggplot2)
ggplot(data_diff[c(1:(49*3),12*49 + 1:(49*3)),],aes(x=month,y=PE,fill=method)) + geom_boxplot()
dev.off()

pdf("Forecast_boxplot1.pdf",width=8,height=4)
library(ggplot2)
ggplot(data_diff,aes(x=month,y=PE,fill=method)) + geom_boxplot()
dev.off()
#MAE_monthly = apply(abs_error, 2, function(x){return(mean(x,na.rm=TRUE))})
# J 1.20, F 1.86, M 1.58, A 0.98, M 1.38, J 1.61, J 1.60, A 1.46
# S 1.14, O 1.43, N 3.33, D 2.03
#mean_monthly = apply(temp_forecast12 - temp_true12, 2, function(x){return(mean(x,na.rm=TRUE))})
# J 0.74, F 1.20, M 1.26, A 0.30, M 0.89, J -0.41, J 0.25, A -0.44
# S 0.10, O -1.08, N 3.305, D -0.20


# 
# Monthchar = c("January","February","March","April","May","June","July", "August", "September", "October", "November", "December")
# #Monthchar = c("July", "August", "September", "October", "November", "December")
# pdf("temp_forecast_error.pdf",width=8,height=5)
# par(mfrow=c(6,2),mar=c(1,2,2,1),oma = c(1, 0, 2, 0))
# for(t in 1:n_forecast){
#     plot((temp_forecast12[,t] - temp_true12[,t])/temp_true12[,t],ylim=c(-1,1), ylab= paste(Monthchar[t],sep=""),
#     type="p", col="red", xlab="")
#     mtext(paste(Monthchar[t],sep=""),side=3,cex=0.5)
# }
# mtext("Relative Errors of Temperature Forecasting", side=3, line=0,outer=T)
# #mtext("Weather Stations", side=1, line=1.5, outer=T)
# dev.off()





# 
# 
# alpha_in12 = parameters$alpha
# char = c("1st", "2nd", "3rd")
# pdf("alpha_forecast.pdf")
# par(mfrow=c(J,1),mar=c(3,1,2,1),oma=c(1,1,1,1))
# 
# for(j in 1:J){
#     K_forecast = matrix(c(parameters$K[1],1,parameters$K[2],0),ncol=2)
#     
#     var.temp = matrix(c(parameters$HJ[j,j],0,0,0),ncol=2)
#     
#     var.i = var.temp
#     var.forecast = var.temp[1,1]
#     for(i in 2:12){
#         var.i = var.temp + K_forecast %*% var.i %*% t(K_forecast)
#         var.forecast = c(var.forecast,var.i[1,1])
#     }
#     plot(1,type="n",xlim=c(0,1200),ylim=c(-50,50),main=paste("Time series of", char[j], "FPC scores",sep=" "))
#     mtext = (side=3)
#     LL = c(alpha_in12[length(alpha_in12)-J+j],alpha_T12[j,] - 1.96 * sqrt(var.forecast))
#     UL = c(alpha_in12[length(alpha_in12)-J+j],alpha_T12[j,] + 1.96 * sqrt(var.forecast))
#     polygon(c(1188:1200,rev(1188:1200)),c(LL,rev(UL)),col="grey75",border=FALSE)
#     lines(c(alpha_in12[seq(j,length(alpha_in12),by=J)],alpha_T12[j,]))
# }
# dev.off()



