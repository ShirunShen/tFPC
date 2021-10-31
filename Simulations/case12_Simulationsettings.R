##set.seed(1)  ##for reproduce purpose.

## domain, triangulation, and data

##mean function mu(x,y)
mufunc <- function(x,y,t){
    #mu <- exp(sqrt(0.5*(x-0.6)^2 + 0.8*y)) * abs(cos(1+t*pi/50));
    mu1 <- 5*exp(sqrt(0.1*x^2 + 0.2*y)) + 5*exp(-sqrt(0.1*x^2 + 0.2*y))
    mu2 <- cos(2*pi*t/12) + 1*t/500 #default 0.05*t/500
    #mu2 <- 1
    return(mu1*mu2)
}

##first PC function qphi1(x,y)
phi1 <- function(x,y){
    phi <- 0.8578*sin(x^2 + 0.5*y^2);
    return(phi)
}

##second PC function phi2(x,y)
phi2 <- function(x,y){
    phi <- 0.8721*sin(0.3*x^2 + 0.6*y^2) - 0.2988*sin(x^2 + 0.5*y^2)
    return(phi)
}

##surface value fi##
fvalue <- function(x,y,i,alpha1,alpha2){
    f <- mufunc(x,y,i)+phi1(x,y)*alpha1[i] + phi2(x,y)*alpha2[i]
    return(f)
}

stations <- function(n,ni){
    ######randomly generated the data points at different time steps 1:n
    ###input:
    ###     n:    number of time steps
    ###     ni:   vector, number of points at each time step
    ###output:
    ###     sur:  a sum(ni) * 2 matrix, collection of randomly generated data points
    sur <- NULL
    for(i in 1:n){
        count = 0
        sur.i = NULL
        while(count < ni[i]){
            xloc <- runif(1,0,2)
            yloc <- runif(1,0,2)
            if(xloc < 0.5 | xloc > 1.5 | yloc < 0.5 | yloc > 1.5){
                sur.i <- rbind(sur.i,c(xloc,yloc))
                count = count + 1
            }
        }
        sur <- rbind(sur, sur.i)
    }
    return(sur)
}

ar2 <- function(a1, a2, sigma, n, n0=200){
    ######AR(2) time series
    ###input:
    ###    a1, a2:   coefficients of time series
    ###    sigma:    variance of Gaussian innovations
    ###    n:        total number of time steps
    ###    n0:       number of burn-in data points
    ###output:
    ###    y:        AR(2) time series
    m = (n0 + n)
    y = rep(0, m)
    w = rnorm(m,0,sigma)
    for(i in 3:m){
        y[i] = a1 * y[i-1] + a2 * y[i-2] + w[i]
    }
    y = y[-c(1:n0)]
    return(y)
}


######useful functions#####

indfunc = function(i.if, n.if){
    ### seize the data from the concreted matrix
    if(i.if==1)  ind=1:n.if[1]
    else         ind=sum(n.if[1:(i.if-1)])+1:n.if[i.if]
    return(ind);
}


