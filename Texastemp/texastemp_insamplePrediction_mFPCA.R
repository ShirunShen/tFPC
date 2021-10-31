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


#station = Texasstation
#boundary = Texasboundary

#pdf("data_texas.pdf",height=4,width=4)
#par(mar=c(4,4,2,2))
#plot(boundary[[1]],type="l",xlab="",ylab="",xlim=c(-107,-93),cex.lab=1.3)
#mtext("Longitude", side=1, line=2)
#mtext("Latitude", side=2, line=2)
#points(station[,1],station[,2],pch=2,cex=0.5)
##plot.t(V,Tri,1,"l")
#dev.off()

#pdf("data_texas.pdf",height=8,width=8)
#par(mfrow=c(2,2),mar=c(0,0,2,0))
#for(i in 1:12){
#    for(j in 1:4){
#plot(boundary[[1]],type="l",xlab="",ylab="",xlim=c(-107,-93),cex.lab=1.3,,main=(j-1)*12+i)
#        mtext("Longitude", side=1, line=2)
#        mtext("Latitude", side=2, line=2)
#        points(station[,1],station[,2],pch=1,cex=.05)
#        text(station[,1],station[,2],labels=round(Texastemp[(j-1)*12+i,],2),cex=.5,col="red")

#    }
#}
##plot.t(V,Tri,1,"l")
#dev.off()


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




#penalty matrix
Ene=energy(V,Tri,d,1)
P=t(omat)%*%Ene%*%(omat)  #with orthogonalizations

J=3    # number of principal component vectors



### mFPCA ###


set.seed(1)
zhat38_mFPCA_summary = NULL
PPE_mFPCA_summary = NULL
MSE_mFPCA = NULL
MAE_mFPCA = NULL


for(ind in 1:12){
    print(ind)
    singlemonth = seq(ind,nrow(Texastemp),by=12)
    Texastemp1 = Texastemp[singlemonth,] #single month
    
    
    ind_singlemonth = NULL
    
    if(singlemonth[1] == 1){
        ind_singlemonth = 1:ni[1]
    }else{
        ind_singlemonth = sum(ni[1:(singlemonth[1]-1)])+1:ni[singlemonth[1]]
    }
    
    
    #pdf("rawdata.pdf")
    #image(as.matrix(Texastemp1))
    #dev.off()
    
    
    ind_station = NULL
    
    z = NULL
    station = NULL
    ni = NULL
    for(i in 1:nrow(Texastemp1)){
        count = 0
        for(j in 1:ncol(Texastemp1)){
            if(!is.na(Texastemp1[i,j])){
                z = c(z,Texastemp1[i,j])
                count = count + 1
                station = rbind(station,Texasstation_38[j,])
                
                ind_station = c(ind_station,j)
            }
        }
        ni = c(ni,count)
    } # need some time!!!
    
    z_mFPCA = z
    
    
    basis = NULL
    for(i in 1:nrow(Texastemp1)){
        for(j in 1:ncol(Texastemp1)){
            if(!is.na(Texastemp1[i,j])){
                basis = rbind(basis,rawbasis2[j,])
            }
        }
    }
    #if(basis[1,1] < 0){basis = -basis}
    
    
    
    
    source("../functions/mFPC/EMinitial_mFPCA_realdata.R")
    Rcpp::sourceCpp("../functions/mFPC/EMalg.cpp")
    
    source("../functions/mFPC/CV_mFPCA_realdata.R")
    source("../functions/mFPC/mFPCA_simplex_realdata.R")
    
    grid_mu = c(10,100,1000,10000)
    grid_pc = c(10,100,1000,10000)
    
    grid_mupc = expand.grid(grid_mu,grid_pc)
    result_penalty=NULL
    for(i in 1:nrow(grid_mupc)){
        result_penalty = c(result_penalty,CrossVal(K_fold=5,J,z,basis,P,ni,grid_mupc[i,1],grid_mupc[i,2]))
    }
    init_mFPCA = as.numeric(unlist(grid_mupc[which.min(result_penalty),]))
    print("the initial penalty parameters ")
    print(init_mFPCA)
    
    
    #init_mFPCA = c(10,10) #default 10 10
    penaltyPar = simplex(K_fold=5,J,z,basis,P,ni,init_mFPCA)
    
    #source("../functions/modelbasedPCA/CV_mFPCA_realdata.R")
    #CrossVal(K_fold=5,J,z,basis,P,ni,1,1)
    
    #penaltyPar = list(c(2,1))
    
    lambmu = penaltyPar[[1]][1]
    lambpc = penaltyPar[[1]][2]
    
    ####
    EMinit = EMinitial(z,ni,basis,J)
    sigma20 = EMinit$sigma20
    D0 = EMinit$D0
    thetamu0 = EMinit$thetamu0
    Theta0 = EMinit$Theta0
    
    EM = EM_algorithm(z,basis,P,ni,lambmu,lambpc,sigma20,D0,thetamu0,Theta0)
    
    thetamu_mFPCA = EM$thetamu
    Thetahat_mFPCA = EM$BigTheta
    
    muhat38_mFPCA = rawbasis238 %*% thetamu_mFPCA
    muhat38_mFPCA = rep(muhat38_mFPCA,length(ni))
    
    temp = rep(0,length(muhat38_mFPCA))
    for(j in 1:J){
        temp = temp + EM$alpha[j,] * as.numeric(rawbasis238 %*% Thetahat_mFPCA[,j])
    }
    zhat38_mFPCA = muhat38_mFPCA + temp
    zhat38_mFPCA_summary = rbind(zhat38_mFPCA_summary,zhat38_mFPCA)
    Pred_error_mFPCA = (zhat38_mFPCA - as.vector(Texastemp38[singlemonth]))
    PPE_mFPCA = (zhat38_mFPCA - as.vector(Texastemp38[singlemonth]))/as.vector(Texastemp38[singlemonth])
    PPE_mFPCA_summary = rbind(PPE_mFPCA_summary,PPE_mFPCA)
    MSE_mFPCA = c(MSE_mFPCA, mean(Pred_error_mFPCA^2,na.rm=TRUE))
   
    MAE_mFPCA = c(MAE_mFPCA, mean(abs(Pred_error_mFPCA),na.rm=TRUE))
}

save.image(file="z38_mFPCA.RData")



write.table(MSE_mFPCA, file="z38mse_mFPCA.txt")

#tMSE_mFPCA = mean((Texastemp38-as.vector(zhat38_mFPCA_summary))^2,na.rm=TRUE)

#write.table(tMSE_mFPCA,file="z38mse_mFPCA.txt",append=TRUE)

write.table(MAE_mFPCA, file="z38mae_mFPCA.txt")

write.table(zhat38_mFPCA_summary,)







# setwd("~/Desktop/tFPC_final/Texastemp/")
# 
# library(Matrix)
# library(MASS)
# library(mgcv)
# library(ape)
# library(fields)
# library(fda)
# load(file="inSample_prediction.RData")



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

