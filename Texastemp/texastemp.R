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



# if change the region of Texas, from here.
#Texasboundary=list(read.table("./Data/texasboundary.dat",header=T))
Texasboundary=list(read.table("./Data/texasregion.txt",header=T))
names(Texasboundary[[1]])=c("x","y")

ind_texas = which(inSide(Texasboundary,x,y))


#pdf("texasmap.pdf")
#plot(Texasboundary[[1]],type="l")
#points(stationlocations[ind_texas,3:2],col="red")
#dev.off()


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

Texastemp = texastemp1[(37*12+1):nrow(texastemp1),] #start from 1935 to 2014, total 80years # 57 means from 1935, #37 means from 1915


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

#            ind_station = c(ind_station,j)
#        }
#    }
#    ni = c(ni,count)
#} # need some time!!!


#write.table(z,file="./Data/temp49_100years.txt")
#write.table(ni,file="./Data/ni49_100years.txt")
#write.table(station,file="./Data/station49_100years.txt")

z = as.matrix(read.table("./Data/temp49_100years.txt"))
ni = as.matrix(read.table("./Data/ni49_100years.txt"))
station = as.matrix(read.table("./Data/station49_100years.txt"))

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


#basis = NULL
#for(i in 1:nrow(Texastemp)){
#    for(j in 1:ncol(Texastemp)){
#        if(!is.na(Texastemp[i,j])){
#            basis = rbind(basis,rawbasis2[j,])
#        }
#    }
#}
#write.table(basis,file="./Data/basis49_100years.txt")

basis = as.matrix(read.table("./Data/basis49_100years.txt"))


## normalize basis
basis1 = B1 %*% omat  #B1 = beval(V,Tri,d,x[ind],y[ind]) # for every time point
nnp = nrow(basis1)

normal.term_spatial = sqrt((t(basis1) %*% basis1)[1,1])

basis = basis/normal.term_spatial
basis1 = basis1/normal.term_spatial





#penalty matrix
Ene=energy(V,Tri,d,1)
P=t(omat)%*%Ene%*%(omat)  #with orthogonalizations



####### tFPCA #######

set.seed(1)
J=3   # number of principal component vectors

splinebasis = create.fourier.basis(c(0,length(ni)), nbasis=51, period = 12)
basismatrix = eval.basis(0:(length(ni)-1), splinebasis)

polybasis = bs(1:length(ni),#knots=c(66)*12,
        degree=3, intercept=FALSE)
basismatrix = cbind(basismatrix,polybasis)


names(basismatrix) = NULL
Pt = fourierpen(splinebasis)

Pt = cbind(Pt, matrix(0,ncol=ncol(basismatrix)-51,nrow=nrow(Pt)))
Pt = rbind(Pt, matrix(0,nrow=ncol(basismatrix)-51,ncol=ncol(Pt)))




nb = ncol(basis)
nc = ncol(basismatrix)
I_nb = diag(nb)
p=2


#### select p ####
#source("../functions/tFPC/AIC.R")
#set.seed(3)
#AIC_tFPCA = NULL
#for(p in 1:3){
#    temp = AIC(p,K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,ni)
#    AIC_tFPCA = c(AIC_tFPCA, temp)
#}

## AIC_tFPCA(1,2,3) =(-6575.928, 140.4404, 69865)
#p = which(AIC_tFPCA == min(AIC_tFPCA))



### test, Oct 26, 2021

source("../functions/tFPC/CV_tFPCA_realdata.R")
grid_mus = c(1,100,10000)
grid_mut = c(0.0001,0.01, 0.1)
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

init = penalty_result[which.min(penalty_result[,1]),-1]


### end test, Oct 26, 2021

set.seed(2)
#init = c(10000,0.0001,10000) #default (10000,0.0001,10000) 

source("../functions/tFPC/tFPCA_simplex_realdata.R")

penalty = simplex(K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,p,ni,init,"EMalgorithm")
#save.image(file="1.RData")


#penalty = list(c(6180, 0.0000618, 6180),1)

lambmus = penalty[[1]][1]
lambmut = penalty[[1]][2]
lambpc  = penalty[[1]][3]

#set.seed(1)
#source("../functions/tFPC/CV_tFPCA_realdata.R")

#CrossVal(K_fold=5,basis,P,basismatrix,Pt,z,V,Tri,length(ni),d,r,J,p,ni,lambmus,lambmut,lambpc,"EMalgorithm_tFPCA0_manifoldthetab") # to check the penalty parameters are reasonable by comparing different parameters setup.


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




system.time({parameters = EMalgorithm(z,basis,basismatrix,P,Pt,ni,HJ0,thetab0,thetac0,
    Theta0,K0,sigma20,lambmus,lambmut,lambpc)})

thetabhat = parameters$thetab
thetachat = parameters$thetac
Thetahat = parameters$Theta


# to present the number of PCs
# pdf("scree.pdf", height=2,width=4)
# par(mar=c(3.5,3.5,3,3))
# plot(diag(parameters$HJ),type="b",pch=19,xlab="",ylab="",ylim=c(0,parameters$HJ[1,1] + 1000),main="Scree Plot")
# mtext("index of pc", side=1, line=2)
# mtext("variance",side=2, line=2)
# dev.off()
# 
# pdf("scree2.pdf",height=4,width=4)
# par(mar=c(3.5,3.5,3,3))
# varHJ = diag(parameters$HJ)
# cumvarHJ = c(varHJ[1],sum(varHJ[1:2]),sum(varHJ[1:3]),sum(varHJ[1:4]))
# plot(cumvarHJ/cumvarHJ[4],type="b",pch=19, xlab="",ylab="",main="Screen Plot",ylim=c(0.5,1.0))
# abline(h=0.99,col="red")
# mtext("index of pc", side=1, line=2)
# mtext("cummulative var %",side=2, line=2)
# dev.off()

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

ym = rep(NA,m*n)

zhat1_tFPCA = zhat1
summary(zhat1_tFPCA)

save.image(file="100years.RData")

# 
# 
# 
# #### plot ######
rm(list=ls())
setwd("~/Desktop/tFPC_final/Texastemp/")

library(Matrix)
library(MASS)
library(mgcv)
library(ape)
library(fields)
library(fda)

load(file="100years.RData")



yearly.mean = NULL
for(i in 38:(nrow(texastemp1)/12)){
    yearly.mean = rbind(yearly.mean, apply(texastemp1[(i-1)*12+1:12,],2,function(x){return(mean(x,na.rm=TRUE))}))
}


mean1_func = function(x){
    temps=0.0
    for(i in 1:length(x)){
        if(is.na(x[i])){temps = NA; break;}
        else{temps = temps + x[i]}
    }
    return(temps/length(x))
}


#### raw data processing ###

## trends
yearly.average = matrix(0, ncol=49,nrow=100)
for(t in 1:100){
    temporary = apply(Texastemp[(t-1)*12+1:12,],2,mean1_func)
    yearly.average[t,] = temporary
}

test = yearly.average
for(j in 1:49){
    test[,j] = yearly.average[,j] - mean(yearly.average[,j],na.rm=TRUE)
}


#### test, Jun 9, 2021
yearly.temp_est = mean(zhat[1:sum(ni[1:12])])
for(i in 2:100){
    indi = sum(ni[1:((i-1)*12)])+ 1:sum(ni[(i-1)*12+1:12])
    yearly.temp_est = c(yearly.temp_est,mean(zhat[indi]))
}


Texastemp_est = matrix(NA,ncol=1200,nrow=49)
ind_na = which(is.na(t(Texastemp)))
Texastemp_est[-ind_na] = zhat
Texastemp_est = t(Texastemp_est)

periodicity.average = matrix(0, ncol=49,nrow=12)
periodicity_est = matrix(0,ncol=49,nrow=12)

for(t in 1:12){
    temporary = apply(Texastemp[t+12*(0:99),],2,function(x){return(mean(x,na.rm=TRUE))})
    periodicity.average[t,] = temporary
    periodicity_est[t,] = apply(Texastemp_est[t+12*(0:99),],2,function(x){return(mean(x,na.rm=TRUE))})
}


temp_spatial = rep(0,nnp)
for(i in 1:length(ni)){
    temp_spatial = temp_spatial + zhat1_tFPCA[(i-1)*nnp + 1:nnp]
}
temp_spatial = temp_spatial/length(ni)



# test
temp_Jan = rep(0,nnp)
temp_Aug = rep(0,nnp)
for(i in seq(1,length(ni),by=12)){
    temp_Jan = temp_Jan + zhat1_tFPCA[(i-1)*nnp+1:nnp]
    temp_Aug = temp_Aug + zhat1_tFPCA[(i+7-1)*nnp+1:nnp]
}
temp_Jan = temp_Jan/(length(seq(1,length(ni),by=12)))
temp_Aug = temp_Aug/(length(seq(1,length(ni),by=12)))


pdf("tFPC_temp.pdf",height=4*2,width=4*2+1)
par(mfrow=c(2,2),mar=c(4,3,2,5),oma=c(2,3,2,2))
plot(1, type='n',main="Yearly-Average Temperature",xlim=c(0,101),ylim=c(10,25),ylab="", xlab="", xaxt="n",yaxt="n", cex.main=1.5)
axis(side=1, at=c(6,26,46,66,86),labels=c(1920,1940,1960,1980,2000),cex.axis=1.2)
axis(side=2,at=c(5,10,15,20,25,30),labels=c(5,10,15,20,25,30),cex.axis=1.2)
mtext("Year",side=1,line=3.0,cex=1.2)
mtext("Temperature in Degree Celsius",side=2,line=3.0,cex=1.2)
for(i in 1:49){
    points(yearly.mean[,i],col=i,cex=0.01)
}
lines(apply(yearly.mean,1,function(x){return(mean(x,na.rm=TRUE))}),col="blue",lwd=2.5,lty=1)
#abline(h=c(20,18,16))
lines(yearly.temp_est,col="red",lwd=2,lty=2)
legend(x=0,y=25,legend=c("tFPC","raw data"),col=c("red","blue"),lty=c(2,1))

#2
periodicity_est.average = apply(periodicity_est,1,mean)

plot(1,type="n",main="Period of Temperature",xaxt="n",#yaxt="n",
     xlab="",ylab="",cex.axis=1.2, cex.main=1.5,xlim=c(1,12),ylim=c(0,35))
mtext("Month",side=1,line=3.0,cex=1.2)
mtext("Temperature in Degree Celsius",side=2,line=3.0,cex=1.2)
for(i in 1:49){
    lines(periodicity.average[,i],col=rgb(255,0,0,80,maxColorValue=255),lty=2)
}
lines(periodicity_est.average,lwd=2.5)
axis(side=1, at=c(1,3,5,7,9,11),labels=c("Jan","Mar","May","Jul","Sep","Nov"),cex.axis=1.2)

#3
fig = matrix(NA,ncol=m,nrow=n)
fig[ind] = temp_Jan
image(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("January Temperature Average", sep=" "),
           zlim=c(-10,50),xaxt="n",cex.axis=1.2,cex.main=1.5)
axis(side=1,at = c(-106,-102,-98,-94),labels=c(-106,-102,-98,-94),cex.axis=1.2)
mtext("Longitude",side=1,line=3.0,cex=1.2)
mtext("Latitude",side=2,line=3.0,cex=1.2)
contour(xm,yn,fig,levels= c(2,6,10)#round(seq(min(fig[ind]),max(fig[ind]),length=5),2)
        , labcex=1.2, add=TRUE)

#4
fig = matrix(NA,ncol=m,nrow=n)
fig[ind] = temp_Aug
image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("August Temperature Average", sep=" "),
           zlim=c(-10,50),xaxt="n",cex.axis=1.2,cex.main=1.5)
axis(side=1,at = c(-106,-102,-98,-94),labels=c(-106,-102,-98,-94),cex.axis=1.2)
mtext("Longitude",side=1,line=3.0,cex=1.2)
mtext("Latitude",side=2,line=3.0,cex=1.2)
contour(xm,yn,fig,levels=c(24,30)
        , labcex=1.2, add=TRUE)
dev.off()





pdf("PrincipalComponent_wholeyear.pdf",height=4, width=4*J+0.1)
par(mfrow=c(1,J),mar=c(4,3,2,5.5),oma=c(2,3,2,2))

ym[ind] = basis1 %*% parameters$Theta[,1]
image.plot(xm,yn,matrix(ym,nrow=m),col=my.colors(50),xlab="",ylab="",main="First PC function",zlim=c(-0.07,0.07),cex.main=1.5,xaxt="n",yaxt="n")
axis(side=1,at=c(-106,-102,-98,-94),labels=c(-106,-102,-98,-94),cex.axis=1.2)
axis(side=2,at=c(26,28,30,32,34,36),labels=c(26,28,30,32,34,36),cex.axis=1.2)
contour(xm,yn,matrix(ym,m,n), levels= round(seq(min(ym[ind]),max(ym[ind]),length=4),2), labcex=1.2, add = TRUE) #round(seq(min(ym[ind]),max(ym[ind]),length=15),2),add=TRUE)

ym[ind] = basis1 %*% parameters$Theta[,2]
image.plot(xm,yn,matrix(ym,nrow=m),col=my.colors(50),xlab="",ylab="",main="Second PC function",zlim=c(-0.07,0.07),cex.main=1.5,xaxt="n",yaxt="n")
axis(side=1,at=c(-106,-102,-98,-94),labels=c(-106,-102,-98,-94),cex.axis=1.2)
axis(side=2,at=c(26,28,30,32,34,36),labels=c(26,28,30,32,34,36),cex.axis=1.2)
contour(xm,yn,matrix(ym,m,n), labcex=1.2, add=TRUE)#levels= c(-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25), labcex=0.8, add= TRUE)
#round(seq(min(ym[ind]),max(ym[ind]),length=15),2),add=TRUE)

ym[ind] = basis1 %*% parameters$Theta[,3]
image.plot(xm,yn,matrix(ym,nrow=m),col=my.colors(50),xlab="",ylab="",main="Third PC function",zlim=c(-0.07,0.07),cex.main=1.5,xaxt="n",yaxt="n")
axis(side=1,at=c(-106,-102,-98,-94),labels=c(-106,-102,-98,-94),cex.axis=1.2)
axis(side=2,at=c(26,28,30,32,34,36),labels=c(26,28,30,32,34,36),cex.axis=1.2)
contour(xm,yn,matrix(ym,m,n), labcex=1.2, add=TRUE)#levels= c(-0.25,-0.2,-0.15,-0.1,-0.05,0.0, 0.05, 0.1,0.15), labcex=0.8, add=TRUE)
#round(seq(min(ym[ind]),max(ym[ind]),length=10),2),add=TRUE)

dev.off()



# pdf("test1.pdf",height=4*6,width=4*2+1)
# par(mfrow=c(6,2),mar=c(4,3,2,5),oma=c(2,3,2,2))
# for(j in 1:12){
#     temp_test = rep(0,nnp)
#     for(i in seq(j,length(ni),by=12)){
#         temp_test = temp_test + zhat1_tFPCA[(i-1)*nnp+1:nnp]
#     }
#     temp_test = temp_test/(length(seq(1,length(ni),by=12)))
#     fig = matrix(NA,ncol=m,nrow=n)
#     fig[ind] = temp_test
#     image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("Estimated Temperature Average", sep=" "),
#                zlim=c(-10,50),xaxt="n",cex.axis=1.2,cex.main=1.5)
#     axis(side=1,at = c(-106,-102,-98,-94),labels=c(-106,-102,-98,-94),cex.axis=1.2)
#     mtext("Longitude",side=1,line=3.0,cex=1.2)
#     mtext("Latitude",side=2,line=3.0,cex=1.2)
# }
# dev.off()
# 
# pdf("texastemp_tFPC.pdf",height=4,width=4*3-0.2)
# par(mfrow=c(1,3),mar=c(4,3,2,3),oma=c(2,3,2,2))
# #1
# plot(1, type='n',main="Estimated Yearly Average Temperature",xlim=c(0,101),ylim=c(5,30),ylab="", xlab="", xaxt="n",yaxt="n", cex.main=1.5)
# axis(side=1, at=c(6,26,46,66,86),labels=c(1920,1940,1960,1980,2000),cex.axis=1.2)
# axis(side=2,at=c(5,10,15,20,25,30),labels=c(5,10,15,20,25,30),cex.axis=1.2)
# mtext("Year",side=1,line=3.0,cex=1.2)
# mtext("Temperature in Degree Celsius",side=2,line=3.0,cex=1.2)
# for(i in 1:49){
#     points(yearly.mean[,i],col=i,cex=0.1)
# }
# lines(apply(yearly.mean,1,function(x){return(mean(x,na.rm=TRUE))}),lwd=2)
# abline(h=c(20,18,16))
# lines(yearly.temp_est,col="red",lwd=2)
# legend(x=0,y=30,legend=c("tFPC estimated","raw data"),col=c("red","black"),lty=c(1,1))
# 
# #2
# periodicity_est.average = apply(periodicity_est,1,mean)
# 
# plot(periodicity_est.average,type="l",main="One Period of Temperature",xaxt="n",#yaxt="n",
#      xlab="",ylab="",cex.axis=1.2, cex.main=1.5,lwd=2.5,ylim=c(0,35))
# mtext("Month",side=1,line=3.0,cex=1.2)
# mtext("Temperature in Degree Celsius",side=2,line=3.0,cex=1.2)
# for(i in 1:49){
#     lines(periodicity.average[,i],col=rgb(255,0,0,80,maxColorValue=255),lty=3)
# }
# axis(side=1, at=c(1,3,5,7,9,11),labels=c("Jan","Mar","May","Jul","Sep","Nov"),cex.axis=1.2)
# 
# 
# #3
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = temp_spatial
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("Estimated Temperature Average", sep=" "),
#            zlim=c(-10,50),xaxt="n",cex.axis=1.2,cex.main=1.5)
# axis(side=1,at = c(-106,-102,-98,-94),labels=c(-106,-102,-98,-94),cex.axis=1.2)
# mtext("Longitude",side=1,line=3.0,cex=1.2)
# mtext("Latitude",side=2,line=3.0,cex=1.2)
# #contour(xm,yn,fig,levels= c(16,20)#round(seq(min(fig[ind]),max(fig[ind]),length=10),2)
# #        , labcex=1.2, add=TRUE)
# 
# dev.off()



#save.image(file="100years.RData")

# pdf("mean_part.pdf")
# plot(1, type='n',xlim=c(0,101),ylim=c(0,30))
# for(i in 1:49){
#     points(yearly.mean[,i],col=i,cex=0.1)
# }
# lines(apply(yearly.mean,1,function(x){return(mean(x,na.rm=TRUE))}),lwd=2)
# abline(h=c(20,18,16))
# lines(yearly.temp_est,col="red",lwd=2)
# dev.off()
# 
# pdf("period_part.pdf")
# periodicity_est.average = apply(periodicity_est,1,mean)
# plot(periodicity_est.average,type="l",lwd=2.5,xaxt="n",xlab="",ylab="",ylim=c(0,35))
# mtext("Month",side=1,line=3.0,cex=1.2)
# mtext("Temperature in Degree Celsius",side=2,line=3.0,cex=1.2)
# for(i in 1:49){
#     lines(periodicity.average[,i],col=rgb(255,0,0,80,maxColorValue=255),lty=3)
# }
# axis(side=1, at=c(1,3,5,7,9,11),labels=c("Jan","Mar","May","Jul","Sep","Nov"),cex.axis=1.2)
# dev.off()
# 
# 
# 
# pdf("mean_spatial_part.pdf")
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = temp_spatial
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("Estimated Temperature Average", sep=" "),
#            zlim=c(-10,50),xaxt="n",cex.axis=1.2,cex.main=1.5)
# axis(side=1,at = c(-106,-102,-98,-94),labels=c(-106,-102,-98,-94),cex.axis=1.2)
# mtext("Longitude",side=1,line=3.0,cex=1.2)
# mtext("Latitude",side=2,line=3.0,cex=1.2)
# contour(xm,yn,fig,levels= round(seq(min(fig[ind]),max(fig[ind]),length=10),2), labcex=1.2, add=TRUE)
# dev.off()

# end test Jun 9, 2021



#### testing, Sep 30, 2020
#temperature.anomalies = as.vector(test)

#polybasis = bs(rep(1:100,49), #seq(10,90,by=10),
#        degree = 3, intercept=FALSE)

#fit <- lm(temperature.anomalies ~ polybasis)

#fit.pred <- polybasis[1:100,1:(ncol(polybasis))] %*% as.vector(fit$coefficients)[1:(ncol(polybasis))]
#fit.pred2 <- predict(fit)

#x = rep(1:100,49)
#x = x[-fit$na.action]

#plot(rep(1:100,49),temperature.anomalies,col="grey",xlab="Year",ylab="temperature anomalies",type="n")
#for(i in 1:49){
#    test1 = test[,i]
#    lines(test1,col=rgb(255,0,0,50,maxColorValue=255),lty=3)
#}
#points(x,fit.pred2,cex=0.2,col="red")
#dev.off()













#test = yearly.average - mean(yearly.average,na.rm=TRUE)

#apply(yearly.average,2,function(x){return(mean(x,na.rm=TRUE))})

#test = yearly.average - apply(yearly.average,2,function(x){return(mean(x,na.rm=TRUE))})

#pdf("test.pdf")
#plot(1,xlim=c(1900,2020),ylim=c(-10,40))
#for(i in 1:49){
#    test1 = ts(test[,i],start=c(1915,1))
#    lines(test1,col=rgb(0,0,255,40,maxColorValue=255))
#    #print(i)
#}
#dev.off()

# spatial mean function
muhat1_spatial = NULL
for(t in 1:length(ni)){
    muhat1_spatial = c(muhat1_spatial, basis1 %*% thetabhat)
}

temp_mean = rep(0, nnp)
for(i in 1:length(ni)){
    temp_mean = temp_mean + muhat1_spatial[(i-1)*nnp + 1:nnp]
}
temp_mean = temp_mean/length(ni)


raw.temp.monthly = apply(Texastemp,1,function(x){return(mean(x,na.rm=TRUE))})

raw.temp.yearly = NULL
for(i in 1:100){
    raw.temp.yearly = c(raw.temp.yearly, mean(raw.temp.monthly[(i-1)*12+1:12]))
}

raw.temp.yearly = raw.temp.yearly - mean(raw.temp.yearly)

## periodicity

periodicity.average = matrix(0, ncol=49,nrow=12)
for(t in 1:12){
    temporary = apply(Texastemp[t+12*(0:99),],2,function(x){return(mean(x,na.rm=TRUE))})
    periodicity.average[t,] = temporary
}

#pdf("test.pdf")
#plot(1,xlim=c(0,13),ylim=c(-10,40))
#for(i in 1:49){
#    lines(periodicity.average[,i],col=i)
#}
#dev.off()

range_period = mean(apply(periodicity.average,2,function(x){return(max(x)-min(x))}))
test2 = periodicity.average - mean(periodicity.average)


pdf("Meanfunction.pdf",height=4,width=4*3-0.2)
par(mfrow=c(1,3),mar=c(4,3,2,3),oma=c(2,3,2,2))

# temp is the estimated trend function
temp = basismatrix[,1] * thetachat[1] + basismatrix[,-(1:51)] %*% thetachat[-(1:51)]#basismatrix[,52] * thetachat[52] + basismatrix[,53] * thetachat[53] + basismatrix[,54] * thetachat[54]
temp = (temp) * mean(muhat1_spatial)
temp = temp - mean(temp)
temp = ts(temp,frequency = 12 , start= c(1915,1))
plot(temp, ylim=c(-5,5),
    type="l",main="Cubic Trend of Mean Function",ylab="", xlab="", xaxt="n",yaxt="n", cex.main=1.5)
axis(side=1, at=c(1920,1940,1960,1980,2000),labels=c(1920,1940,1960,1980,2000),cex.axis=1.2)
axis(side=2,at=c(-4,-2,0,2,4),labels=c(-4,-2,0,2,4),cex.axis=1.2)
mtext("Year",side=1,line=3.0,cex=1.2)
mtext("Temperature in Degree Celsius",side=2,line=3.0,cex=1.2)
#lines(1915:2014,raw.temp.yearly,col="red",lty=3)
legend(1915,5,legend=c("Yearly Average of Raw Data","tFPCA Trend Function"),col=c("red","black"),lty=c(3,1),cex=1.2)
for(i in 1:49){
    test1 = ts(test[,i],start=c(1915,1))
    lines(test1,col=rgb(255,0,0,50,maxColorValue=255),lty=3)
}
lines(temp,type="l",lwd=2.5)


temp0 = -basismatrix[1:12,2:51] %*% thetachat[2:51]
temp1 = temp0/(max(temp0)-min(temp0)) * range_period
plot(temp1,type="l",main="One Period of Mean Function",xaxt="n",#yaxt="n",
xlab="",ylab="",cex.axis=1.2, cex.main=1.5,lwd=2.5,ylim=c(-16,12))
mtext("Month",side=1,line=3.0,cex=1.2)
mtext("Temperature in Degree Celsius",side=2,line=3.0,cex=1.2)
for(i in 1:49){
    lines(test2[,i],col=rgb(255,0,0,80,maxColorValue=255),lty=3)
}
axis(side=1, at=c(1,3,5,7,9,11),labels=c("Jan","Mar","May","Jul","Sep","Nov"),cex.axis=1.2)


fig = matrix(NA,ncol=m,nrow=n)
fig[ind] = -temp_mean
image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("Mean function (spatial effects)", sep=" ") ,zlim=c(-0.03,0.07),xaxt="n",cex.axis=1.2,cex.main=1.5)
axis(side=1,at = c(-106,-102,-98,-94),labels=c(-106,-102,-98,-94),cex.axis=1.2)
mtext("Longitude",side=1,line=3.0,cex=1.2)
mtext("Latitude",side=2,line=3.0,cex=1.2)
contour(xm,yn,fig,levels= round(seq(min(fig[ind]),max(fig[ind]),length=10),2), labcex=1.2, add=TRUE)


dev.off()


#save.image(file="100years.RData")




# temp_mean_t = -basismatrix %*% thetachat #/ mean(basismatrix %*% thetachat)
# pdf("mean_of_time_tFPCA.pdf")
# #pdf("test2.pdf")
# temp_mean_t = ts(temp_mean_t, frequency = 12 , start= c(1915,1))
# plot(temp_mean_t,main="Mean function (temporal effects)", ylab="Mean function mu2(t) without normalization", xlab="Time", cex.lab=1.2, cex.main=1.5)
# dev.off()
# 
# 
# 
# temp_mean = rep(0, nnp)
# for(i in 1:length(ni)){
#     temp_mean = temp_mean + zhat1[(i-1)*nnp + 1:nnp]
# }
# temp_mean = temp_mean/length(ni)
# 
# 
# pdf("average_temperature.pdf")
# #pdf("test1.pdf")
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = temp_mean
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("Estimated Temperature Function", sep=" ") ,zlim=c(-10,40),cex.main=1.5)
# contour(xm,yn,fig, levels= c(-5,0,5,10,15,20,25,30,35),add=TRUE) #round(seq(min(fig[ind]),max(fig[ind]),length=10),2), labcex=1.2, add=TRUE)
# dev.off()
# 
# 
# 
# z_Jan = rep(0,nnp)
# z_Apr = rep(0,nnp)
# z_Jul = rep(0,nnp)
# z_Oct = rep(0,nnp)
# 
# n_mon = nrow(basismatrix)/12
# for(i in 1:n_mon){
#     z_Jan = z_Jan + muhat1[(i-1)*12*nnp+0*nnp+1:nnp]
#     z_Apr = z_Apr + muhat1[(i-1)*12*nnp+3*nnp+1:nnp]
#     z_Jul = z_Jul + muhat1[(i-1)*12*nnp+6*nnp+1:nnp]
#     z_Oct = z_Oct + muhat1[(i-1)*12*nnp+9*nnp+1:nnp]
# }
# z_Jan = z_Jan/n_mon
# z_Apr = z_Apr/n_mon
# z_Jul = z_Jul/n_mon
# z_Oct = z_Oct/n_mon
# 
# 
# pdf("average_mean_function.pdf",height=4, width=4*4+2)
# par(mfrow=c(1,4),mar=c(2,3,2,5.5))
# 
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = z_Jan
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("January Mean Function", sep=" ") ,zlim=c(-10,40),cex.main=1.5)
# #contour(xm,yn,fig, levels= c(5,8,10),add=TRUE)
# 
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = z_Apr
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("April Mean Function", sep=" ") ,zlim=c(-10,40),cex.main=1.5)
# #contour(xm,yn,fig, levels= c(10,12,15,18),add=TRUE)
# 
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = z_Jul
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("July Mean Function", sep=" ") ,zlim=c(-10,40),cex.main=1.5)
# #contour(xm,yn,fig, levels= c(10,15,20,25),add=TRUE)
# 
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = z_Oct
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("October Mean Function", sep=" ") ,zlim=c(-10,40),cex.main=1.5)
# #contour(xm,yn,fig, levels= c(5,10,15),add=TRUE)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# z_Jan = rep(0,nnp)
# z_Apr = rep(0,nnp)
# z_Jul = rep(0,nnp)
# z_Oct = rep(0,nnp)
# 
# n_mon = nrow(basismatrix)/12
# for(i in 1:n_mon){
#     z_Jan = z_Jan + zhat1[(i-1)*12*nnp+0*nnp+1:nnp]
#     z_Apr = z_Apr + zhat1[(i-1)*12*nnp+3*nnp+1:nnp]
#     z_Jul = z_Jul + zhat1[(i-1)*12*nnp+6*nnp+1:nnp]
#     z_Oct = z_Oct + zhat1[(i-1)*12*nnp+9*nnp+1:nnp]
# }
# z_Jan = z_Jan/n_mon
# z_Apr = z_Apr/n_mon
# z_Jul = z_Jul/n_mon
# z_Oct = z_Oct/n_mon
# 
# 
# pdf("Temperature_function.pdf",height=4, width=4*4+2)
# par(mfrow=c(1,4),mar=c(2,3,2,5.5))
# 
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = z_Jan
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("January Function", sep=" ") ,zlim=c(-10,40),cex.main=1.5)
# contour(xm,yn,fig, levels= c(-5,0,2.5,5,7.5,10,12.5,15),add=TRUE,labcex=1.0)
# 
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = z_Apr
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("April Function", sep=" ") ,zlim=c(-10,40),cex.main=1.5)
# contour(xm,yn,fig, levels= c(10,15,20,22,25,30,35),add=TRUE,labcex=1.0)
# 
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = z_Jul
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("July Function", sep=" ") ,zlim=c(-10,40),cex.main=1.5)
# contour(xm,yn,fig, levels= c(20,25,30,35),add=TRUE,labcex=1.0)
# 
# fig = matrix(NA,ncol=m,nrow=n)
# fig[ind] = z_Oct
# image.plot(xm,yn,fig,col=my.colors(50),xlab="",ylab="",main=paste("October Function", sep=" ") ,zlim=c(-10,40),cex.main=1.5)
# contour(xm,yn,fig, levels= c(0,5,10,15,17.5,20,25,30,35),add=TRUE,labcex=1.0)
# dev.off()
# 
# 
# 
# 
# 
# 
# #library(animation)
# #oopt = ani.options(interval=1)
# #saveGIF({
# #    for(i in 1:length(ni)){
# #        par(mfrow=c(1,1),pty="s")
# #        fig1 = matrix(NA,ncol=m,nrow=n)
# #        fig1[ind] = zhat1[(i-1)*nnp+1:nnp]
# #        image.plot(xm,yn,fig1,xlab="",ylab="",main=paste("estimated mean function", i, sep=""),zlim=c(-5,40))#,col = my.colors(50)   )#,zlim=c(-5,38))
# # contour(xm,yn,fig1,levels=round(seq(min(fig1[ind]),max(fig1[ind]),length=10),2),add=TRUE)
#         
# #        ani.pause()
# #    }
# #},ani.weight=400*1,ani.height=400+20,movie.name="TexasTempPenalty.gif")
# #ani.options(oopt)
# 



