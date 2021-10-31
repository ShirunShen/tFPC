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


#write.table(z,file="./Data/temp.txt")
#write.table(ni,file="./Data/ni.txt")
#write.table(station,file="./Data/station.txt")

z = as.matrix(read.table("./Data/temp49_100years.txt"))

ni = as.matrix(read.table("./Data/ni49_100years.txt"))
station = as.matrix(read.table("./Data/station49_100years.txt"))
 
z_mFPCA = z



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



#penalty matrix
Ene=energy(V,Tri,d,1)
P=t(omat)%*%Ene%*%(omat)  #with orthogonalizations

J=3


#### mFPCA for each single month ####

set.seed(1)
temp_forecast = NULL
for(ind in 1:n_forecast){
    print(ind)
    singlemonth = seq(ind,nrow(Texastemp),by=12)
    Texastemp1 = Texastemp[singlemonth,]
    

    
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
                station = rbind(station,Texasstation[j,])
                
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
    
    
    temp_forecast = cbind(temp_forecast, rawbasis2 %*% thetamu_mFPCA)
}


temp_true12 = t(Texastemp_out)

#write.table(temp_true12, file="true2014.txt")
write.table(temp_forecast,file="mFPCA2014.txt")

save.image("forecast_mFPCA.RData")
    



