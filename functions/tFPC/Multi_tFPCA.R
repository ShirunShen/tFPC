my.colors=function (n, alpha = 1)
{
    if ((n <- as.integer(n[1L])) > 0) {
        j <- n%/%4
        i <-  n-j
        rev(c(rainbow(i, start = 0, end = 1/6, alpha = alpha), if (j >
        0) hsv(h = 1/6, s = seq.int(from = 1 - 1/(2 * j),
        to = 1/(2 * j),
        length.out = j), v = 1, alpha = alpha)))
    }
    else character(0L)
}


indfunc = function(i.if, n.if){
  ### seize the data from the concreted matrix
  if(i.if==1)  ind=1:n.if[1]
  else         ind=sum(n.if[1:(i.if-1)])+1:n.if[i.if]
  return(ind);
}


emp.area=function(fsb,n.area){
### calculate the area of the domain using Monte Carlo method
    lowlim=100*.Machine$double.eps #define this in the function
    set.seed(20130425) #use fixed seed

    fsb.x.lb=min(fsb[[1]]$x)-lowlim
    fsb.x.ub=max(fsb[[1]]$x)+lowlim
    fsb.y.lb=min(fsb[[1]]$y)-lowlim
    fsb.y.ub=max(fsb[[1]]$y)+lowlim
    x=runif(n.area,fsb.x.lb,fsb.x.ub)
    y=runif(n.area,fsb.y.lb,fsb.y.ub)
    return(sum(inSide(fsb,x,y))/n.area*((fsb.x.ub-fsb.x.lb)*(fsb.y.ub-fsb.y.lb)))
}

emp.sample = function(fsb,nsample,lowlim=100*.Machine$double.eps){
### sample locations within boundary "fsb"

        set.seed(20130425) #use fixed seed

	fsb.x.lb=min(fsb[[1]]$x)-lowlim
	fsb.x.ub=max(fsb[[1]]$x)+lowlim
	fsb.y.lb=min(fsb[[1]]$y)-lowlim
	fsb.y.ub=max(fsb[[1]]$y)+lowlim
	### get sample locations
	i=1
	x={}
	y={}
	while(i<=nsample){
		xx=runif(1,fsb.x.lb,fsb.x.ub)
		yy=runif(1,fsb.y.lb,fsb.y.ub)
		if(inSide(fsb,xx,yy)) {
			x=c(x,xx)
			y=c(y,yy)
			i=i+1
		}
	}
	return(cbind(x,y))
}

emp.inner=function(f1,f2,fsb,area){
### calculate empirical integration of f1*f2 over domain with boundary fsb using Monte Carlo method, 2000 random points
	loc.inner=emp.sample(fsb,2000)
	return( sum(f1(loc.inner[,1],loc.inner[,2])*f2(loc.inner[,1],loc.inner[,2])*area/nrow(loc.inner)))
}



beval=function(v.se,t.se,d.se,x.se,y.se){
    #create b-spline basis
    #v.se: vertices
    #t.se: edges
    #d.se: bspline degree
    #x.se, y.se: (x,y) coordinates
	tol=100*.Machine$double.eps
	ntri.se=nrow(t.se)
	nbasis.se=(d.se+1)*(d.se+2)/2
	ind.se=rep(1,length(x.se))
	bmatrix.se=matrix(0,length(x.se),nbasis.se*nrow(t.se))
#        browser()
	for (i in 1:ntri.se){
      #barycentric coordinates of (x.se, y.se) within triangle i
      result=bary(v.se[t.se[i,1],],v.se[t.se[i,2],],v.se[t.se[i,3],],x.se,y.se)
      lam1.se=result[,1];lam2.se=result[,2];lam3.se=result[,3];
      I.se=((lam1.se>=-tol) *(lam2.se>=-tol) * (lam3.se>=-tol)==1)
      I.se=which(I.se&ind.se) #only points within triangle i and no bspline coordinates
      if(length(I.se)>0){ #(x.se,y.se) have >0 barycentric coordinates
          ind.se[I.se]=0#to avoid double count on-triangle-edge points
		    for(j in 1:nbasis.se){
            c.se=rep(0,nbasis.se)
            c.se[j]=1
            bmatrix.se[I.se,(i-1)*nbasis.se+j]=loceval(lam1.se[I.se],lam2.se[I.se],lam3.se[I.se],c.se)
		    }
      }
  }
    return(bmatrix.se)
}


################### loceval ###################
loceval=function(lam1.le,lam2.le,lam3.le,bcoef.le){
	nc.le=length(bcoef.le)
	d.le=degree(nc.le)
#        browser()
	for (j in 1:d.le) bcoef.le=de.cast.step(lam1.le,lam2.le,lam3.le,bcoef.le)
	y=t(bcoef.le)
	return(y)
	}


################## bary ################## 


bary=function(v1.ba,v2.ba,v3.ba,x.ba,y.ba){
    #return bary centric coordinates
    one.ba=matrix(1,1,length(x.ba))
	A.ba=rbind(rep(1,3),c(v1.ba[1],v2.ba[1],v3.ba[1]),c(v1.ba[2],v2.ba[2],v3.ba[2]))
	lam.ba=solve(A.ba) %*% rbind(one.ba,t(x.ba),t(y.ba))
	return(t(as.matrix(lam.ba)))
	}

####################### de.cast.step
#######################

de.cast.step=function(lam1.bo,lam2.bo,lam3.bo,Bin.bo){
    #generate spline basis
	#lam1.bo=as.matrix(lam1.bo)
	#lam2.bo=as.matrix(lam2.bo)
	#lam3.bo=as.matrix(lam3.bo)
	Bin.bo=as.matrix(Bin.bo)
	m.bo=nrow(Bin.bo)
	d.bo=degree(m.bo)
	n.bo=length(lam1.bo)
	result=indices(d.bo)
	I.bo=result[,1];J.bo=result[,2];K.bo=result[,3]
	result=indices(d.bo-1)
	I1.bo=result[,1];J1.bo=result[,2];K1.bo=result[,3]
	index1.bo=locate(cbind(I1.bo+1,J1.bo,K1.bo),cbind(I.bo,J.bo,K.bo))
	index2.bo=locate(cbind(I1.bo,J1.bo+1,K1.bo),cbind(I.bo,J.bo,K.bo))
	index3.bo=locate(cbind(I1.bo,J1.bo,K1.bo+1),cbind(I.bo,J.bo,K.bo))
	if (ncol(Bin.bo)==1)
            Bout.bo=Bin.bo[index1.bo]%*%t(lam1.bo)+Bin.bo[index2.bo]%*%t(lam2.bo)+Bin.bo[index3.bo]%*%t(lam3.bo)
        else{
            if(length(lam1.bo)>1){
                Bout.bo=Bin.bo[index1.bo,]%*%diag(lam1.bo,length(lam1.bo),length(lam1.bo))+Bin.bo[index2.bo,]%*%diag(lam2.bo,length(lam2.bo),
                length(lam2.bo))+Bin.bo[index3.bo,]%*%diag(lam3.bo,length(lam3.bo),length(lam3.bo))}
            else{
		Bout.bo=Bin.bo[index1.bo,]*lam1.bo+Bin.bo[index2.bo,]*lam2.bo+Bin.bo[index3.bo,]*lam3.bo}
        }
	return(Bout.bo)
	}

##################### degree
#####################

degree=function(m.dg){
    #calculate degree of bsplines
    #m.dg: number of basis (m=(d+1)(d+2)/2)
    d.dg=(-3+sqrt(8*m.dg+1))/2
    return(d.dg)
}

##################### indices
#####################

indices=function(d.in){
    #generate indices i,j,k: i+j+k=d.in
	m.in=(d.in+1)*(d.in+2)/2
	I.in=rep(0,m.in)
	J.in=I.in
	K.in=I.in
	Mark.in=1
	for (j in seq(d.in,0,-1)){
		I.in[Mark.in:(Mark.in+j)]=seq(j,0,-1)
		J.in[Mark.in:(Mark.in+j)]=0:j
		K.in[Mark.in:(Mark.in+j)]=(d.in-j)*rep(1,j+1)
		Mark.in=Mark.in+j+1
		}
		return(cbind(I.in,J.in,K.in))
	}


#################### locate
####################
locate=function(matrix1.lc,matrix2.lc){
	colnames(matrix1.lc)=NULL
	colnames(matrix2.lc)=NULL
	n1.lc=nrow(matrix1.lc)
	n2.lc=nrow(matrix2.lc)
	ind.lc=rep(0,n1.lc)
	for(j in 1:n1.lc){
		for(k in 1:n2.lc){
			if(all(matrix1.lc[j,]==matrix2.lc[k,])) {ind.lc[j]=k;break}
			}
		}
	return(ind.lc)
    }

################### plot.t
###################

 plot.t= function(v.pt,t.pt,option.t,lp){
     #v.pt: vertices
     #t.pt: triangle vertice indices based on v.pt [row1: triangle 1]
 	ntri.pt=nrow(t.pt)
  	tri.pt=matrix(0,ntri.pt,6)
 	for (j in 1:ntri.pt) {
 		tri.pt[j,]=c(v.pt[t.pt[j,1],],v.pt[t.pt[j,2],],v.pt[t.pt[j,3],])
 		}
 	x.pt=t(tri.pt[,c(seq(1,5,2),1)])
 	y.pt=t(tri.pt[,c(seq(2,6,2),2)])
#        browser()
 	if(lp=="l") #add lines to current plot
  	lines(x.pt,y.pt,col="blue",pch=20)
  	if(lp=="p")
   	plot(x.pt,y.pt,col="blue",pch=20)
  	segments(x.pt[1,],y.pt[1,],x.pt[2,],y.pt[2,],col="blue")
  	segments(x.pt[2,],y.pt[2,],x.pt[3,],y.pt[3,],col="blue")
  	segments(x.pt[1,],y.pt[1,],x.pt[3,],y.pt[3,],col="blue")
  	if (option.t==2) {
  		centerX=apply(x.pt[1:3,],2,mean)
  		centerY=apply(y.pt[1:3,],2,mean)
  		text(centerX,centerY,as.character(1:dim(x.pt)[2]))
  		}
 	}

 ################## plot.z
 ##################
 plot.z=function(v.pz,t.pz,z.pz){
 	ntri.pz=dim(t.pz)[1]
 	tri.pz=matrix(0,ntri.pz,9)
    x.pz=t(tri.pz[,c(seq(1,9,3),1)])
    y.pz=t(tri.pz[,c(seq(2,9,3),2)])
    z.pz=t(tri.pz[,c(seq(3,9,3),3)])??
    xmax.pz=max(x.pz);xmin=min(x.pz)
    ymax.pz=max(y.pz);ymin=min(y.pz)
    zmax.pz=max(z.pz);zmin=min(z.pz)
    scatterplot3d(x.pz,y.pz,z.pz,col="blue")
 	}



 ################## tdata
 ################## return the vertex, edge and etc
 tdata=function(v.t,t.t){
     #v.t: vertices coordinates ((x,y))
     #t.t: triangle vertices indices (counter clockwise)
     ntri.t=nrow(t.t) #number of triangles
     #initial empty space
     Edges.t=matrix(rep(integer(0),4),ncol=2) #edge's end points indices in V
     TE.t=rep(integer(0),ntri.t) #indicate triangle's edges; row: triangle; col: edge
     numEdges.t=0

     for (j in 1:ntri.t){
         Tj.t=t.t[j,]
         for (k in 1:3){
             edge.t=c(min(Tj.t[k],Tj.t[k%%3+1]),max(Tj.t[k],Tj.t[k%%3 + 1])) #edge.t=sort(c(Tj.t[k],Tj.t[(k+1)%%3])) #ordered vertices for edge.t
             if (nrow(Edges.t)>0 )
             {edgenum.t=which((edge.t[1]==Edges.t[,1])*(edge.t[2]==Edges.t[,2])==1)}
             else
             {edgenum.t=integer(0)}

             if (length(edgenum.t)==0) {
                 Edges.t=rbind(Edges.t,edge.t)
                 numEdges.t=numEdges.t+1
                 TE.t=cbind(TE.t, rep(0,ntri.t)) #start from a boundary edge
                 edgenum.t=numEdges.t
             }
             TE.t[j,edgenum.t]=1
         }#end k loop
     }#end j loop

     numV.t=dim(v.t)[1]
     TV.t=matrix(0,ntri.t,numV.t) #indicate triangle vertices; row: triangle;col: vertices
     for (j in 1:ntri.t){
         TV.t[j,t.t[j,]]=1
     }

     EV.t=matrix(0,numEdges.t,numV.t) #indicate edge vertices -- row: edge; col: vertixes
     for (j in 1:numEdges.t){
         EV.t[j,Edges.t[j,]]=1
     }
    #bdr.t=findbdt(t.t,v.t,edge.t,te.t,ev.t)

     return(list(Edges.t,TE.t,TV.t,EV.t))
 }


 ################## triarea
 ##################
 triarea=function(v1.ta,v2.ta,v3.ta){
    #calculate area of a triangle defined by the 3 vertices
 	x.ta=v1.ta[1];y.ta=v1.ta[2]
 	a.ta=v2.ta[1];b.ta=v2.ta[2]
 	c.ta=v3.ta[1];d.ta=v3.ta[2]
 	A.ta=(a.ta-x.ta)*(d.ta-y.ta)-(c.ta-x.ta)*(b.ta-y.ta)
 	A.ta=A.ta/2
 	return(A.ta)
 	}

 ################## crarrays
 ##################
 crarrays=function(d.cr,r.cr){
 	I1.cr=list();I2.cr=list()
 	for (j in  0:r.cr){
 		result=crcellarrays(d.cr,j)
       I1.cr[[j+1]]=result[[1]]
       I2.cr[[j+1]]=result[[2]]
 		}
   return(list(I1.cr,I2.cr))
 	}

 ################## flip
 ##################
 flip=function(matrix.fl){
     #flip.mat=matrix.f1[nrow(matrix.f1):1,]
 	matrix.fl=as.matrix(matrix.fl)
 	n.fl=nrow(matrix.fl)
 	m.fl=matrix(0,n.fl,n.fl)
 	for(j in 1:n.fl){
 	m.fl[j,n.fl-j+1]=1}
 	return(m.fl%*%matrix.fl)
 }                           #turn the matrix upside down
 ##################crcellarrays
 ##################
 crcellarrays=function(d.cs,r.cs){
 	I1.cs=list();I2.cs=list();
 	result=cr.indices(d.cs,r.cs)
 	J1.cs=result[[1]];J2.cs=result[[2]]
 	D1.cs=rep(0,d.cs+1)
 	D2.cs=rep(0,d.cs+1)
 	s1.cs=d.cs+1
 	s2.cs=1
 	for (j in 1:(d.cs+1)){
 		D1.cs[j]=s1.cs
 		s1.cs=s1.cs+d.cs+1-j
 		D2.cs[j]=s2.cs
 		s2.cs=s2.cs+d.cs+2-j
 		}
 	I2.cs[[1]]=flip(J2.cs)
 	Temp.cs=D1.cs-r.cs
 	I2.cs[[2]]=flip(Temp.cs[1:(d.cs+1-r.cs)])
 	Temp.cs=D2.cs+r.cs
 	I2.cs[[3]]=as.matrix(Temp.cs[1:(d.cs+1-r.cs)])
 	I1.cs[[1]]=J1.cs
 	Temp.cs=matrix(0,nrow(J1.cs),ncol(J1.cs))
 	for (j in 0:r.cs){
 		Temp.cs[,j+1]=D1.cs[(j+1):(d.cs+1-r.cs+j)]
 		}
 	loc.cs=r.cs+2
 	back.cs=r.cs+1
 	if(r.cs>0){
 	for (j in 1:r.cs){
 		for(k in 0:(r.cs-j)){
 			Temp.cs[,loc.cs]=Temp.cs[,loc.cs-back.cs]-1
 			loc.cs=loc.cs+1
 			}
 		back.cs=back.cs-1
 		}
 		}
 	I1.cs[[2]]=Temp.cs
 	Temp.cs=matrix(0,nrow(J1.cs),ncol(J1.cs))
 	for(j in 0:r.cs){
 		Temp.cs[,j+1]=D2.cs[(j+1):(d.cs+1-r.cs+j)]
 		}
 	loc.cs=r.cs+2
 	back.cs=r.cs+1
 	if(r.cs>0){
 	for (j in 1:r.cs){
 		for(k in 0:(r.cs-j)){
 			Temp.cs[,loc.cs]=Temp.cs[,loc.cs-back.cs]+1
 			loc.cs=loc.cs+1
 			}
 		back.cs=back.cs-1
 		}}
 		I1.cs[[3]]=flip(Temp.cs)
 		return(list(I1.cs,I2.cs))
 	}

 ################## cr.indices
 ##################
 cr.indices=function(d.ci,r.ci){
 	I1.ci=integer(0)
 	start.ci=1
 	D.ci=d.ci+1
 	for (j in 0:r.ci){
 		for (k in 0:(r.ci-j)){
 			new.col.ci=(start.ci+k):(start.ci+k+d.ci-r.ci)
 			I1.ci=cbind(I1.ci,new.col.ci)
 			}
 			start.ci=start.ci+D.ci
 			D.ci=D.ci-1
 		}
    I2.ci=(-r.ci*r.ci/2+r.ci*(d.ci+3/2)+1):(-r.ci*r.ci/2+r.ci*(d.ci+1/2)+d.ci+1)
    return(list(as.matrix(I1.ci),as.matrix(I2.ci)))
  	}

 ################## smoothness
 ##################

 smoothness=function(V.sm,T.sm,d.sm,r.sm){
     #define smoothness matrix H
    result=tdata(V.sm,T.sm)  ##tdata:store triangle in 4 different way
 	E.sm=result[[1]];TE.sm=result[[2]];TV.sm=result[[3]];EV.sm=result[[4]]
 	int.sm=which(apply(TE.sm,2,sum)>1) #interior edge indicator
    n.sm=length(int.sm)                #number of interior edges
    N.sm=0          #number of constraint point of each interior edge
    Neq.sm=0        #number of smoothness constraints on each interior edge
    result=crarrays(d.sm,r.sm)   #need to read again
 	I1.sm=result[[1]]
 	I2.sm=result[[2]]
  	I.sm=list();J.sm=list();K.sm=list()
 	for(j in 0:r.sm){
 		N.sm=N.sm+((j+1)*(j+2)/2+1)*(d.sm+1-j)
 		Neq.sm=Neq.sm+d.sm+1-j
 		result=indices(j)
 		LI.sm=result[,1]
 		LJ.sm=result[,2]
 		LK.sm=result[,3]
 		I.sm[[j+1]]=LI.sm
 		J.sm[[j+1]]=LJ.sm
 		K.sm[[j+1]]=LK.sm
 		}
 	nbasis.sm=(d.sm+1)*(d.sm+2)/2
 	Index1.sm=matrix(0,N.sm*n.sm,1)
 	Index2.sm=matrix(0,N.sm*n.sm,1)
 	values.sm=Index1.sm
 	A.sm=matrix(c(1,2,2,3,3,1,2,1,3,2,1,3),ncol=2,byrow=T)

#        browser()

 	for (j in 1:n.sm){
        k.sm=int.sm[j]                 #int.sm:indice of interior edges
 		AdjT.sm=which(TE.sm[,k.sm]!=0) #triangles with edge k.sm

        v1.sm=E.sm[k.sm,1]             #vertex 1 of the edge k.sm
        t1.sm=AdjT.sm[1]               #triangle 1 adjacent to edge k.sm
        T1.sm=T.sm[AdjT.sm[1],]        #vertexes of the triangle 1 ad...
        v2.sm=E.sm[k.sm,2]             #vertex 2 of the edge k.sm
        t2.sm=AdjT.sm[2]               #triangle 2 adjacent to edge k.sm
        T2.sm=T.sm[AdjT.sm[2],]        #vertexes of the triangle 2 ad...

        i1.sm=which(T1.sm==v1.sm)      #indice of vertex of triangle 1 that start the edge we consider now
        j1.sm=which(T2.sm==v1.sm)      #indice of vertex of triangle 2 that.
        i2.sm=which(T1.sm==v2.sm)      #indice of vertex of triangle 1 that end the edge we consider now
        j2.sm=which(T2.sm==v2.sm)      #indice of vertex of triangle 2 that.

 		e1.sm = locate(matrix(c(i1.sm,i2.sm),nrow=1),A.sm)
 		e2.sm = locate(matrix(c(j1.sm,j2.sm),nrow=1),A.sm)

 		if (e1.sm >3) {
 			e1.sm=e1.sm-3
 			temp.sm=T1.sm
 			T1.sm=T2.sm
 			T2.sm=temp.sm
 			temp.sm=e1.sm
 			e1.sm=e2.sm
 			e2.sm=temp.sm
 			temp.sm=t1.sm
 			t1.sm=t2.sm
 			t2.sm=temp.sm
 			} else
            e2.sm=e2.sm-3              #why -3???

 		v4.sm=T2.sm[T2.sm%in%c(v1.sm,v2.sm)==FALSE]
        V4.sm=V.sm[v4.sm,]             #coordinate of v4.sm
    result=bary(V.sm[T1.sm[1],],V.sm[T1.sm[2],],V.sm[T1.sm[3],],V4.sm[1],V4.sm[2])
 		lam1.sm=result[,1]
 		lam2.sm=result[,2]
 		lam3.sm=result[,3]
 		lambda.sm=c(lam1.sm,lam2.sm,lam3.sm)

 		if(e1.sm==2){
 			temp.sm=lambda.sm[1]
 			lambda.sm[1]=lambda.sm[2]
			lambda.sm[2]=lambda.sm[3]
			lambda.sm[3]=temp.sm
 			}
 		if(e1.sm==3){
 			temp.sm=lambda.sm[2]
 			lambda.sm[2]=lambda.sm[3]
 			lambda.sm[3]=temp.sm
 			}
 		VarCT.sm=0
 		EqCt.sm=0
 		for (k in 0:r.sm){
 			lam.sm=factorial(k)/(gamma(I.sm[[k+1]]+1)*gamma(J.sm[[k+1]]+1)*gamma(K.sm[[k+1]]+1))*lambda.sm[1]^I.sm[[k+1]]*lambda.sm[2]^J.sm[[k+1]]*lambda.sm[3]^K.sm[[k+1]]
 			T1mat.sm=as.matrix(I1.sm[[k+1]][[e1.sm]])
            T2vector.sm=I2.sm[[k+1]][[e2.sm]]  ##problem located: e2.sm=-2!!!
 			numeq.sm=nrow(T1mat.sm)
 			numVar.sm=(ncol(T1mat.sm)+1)*numeq.sm
 			T1Values.sm=matrix(1,nrow(T1mat.sm),ncol(T1mat.sm))%*%diag(lam.sm)
 			eqnums.sm=(j-1)*(Neq.sm)+EqCt.sm+diag(1:numeq.sm)%*%matrix(1,numeq.sm,ncol(T1mat.sm)+1)
 			Index1.sm[((j-1)*N.sm+VarCT.sm+1):((j-1)*N.sm+VarCT.sm+numVar.sm)]=eqnums.sm
 			Index2.sm[((j-1)*N.sm+VarCT.sm+1):((j-1)*N.sm+VarCT.sm+numVar.sm)]=c((t1.sm-1)*nbasis.sm+T1mat.sm,(t2.sm-1)*nbasis.sm+T2vector.sm)
 			values.sm[((j-1)*N.sm+VarCT.sm+1):((j-1)*N.sm+VarCT.sm+numVar.sm)]=c(T1Values.sm,-1*rep(1,length(T2vector.sm)))
 			VarCT.sm=VarCT.sm+numVar.sm
 			EqCt.sm=EqCt.sm+numeq.sm
 			}
 		}
 			H.sm=matrix(0,n.sm*Neq.sm,nrow(T.sm)*nbasis.sm)
 			H.sm[cbind(Index1.sm,Index2.sm)]=values.sm
 			return(H.sm)
 	}

 ###################### energy
 ###################### evaluate energy matrix
 energy=function(V.e,T.e,d.e,index.e){
     #penalty matrix
 	ntri.e=nrow(T.e)
    D.e=(d.e+1)*(d.e+2)/2    #number of i+j+k = d
 	Dsq.e=D.e^2
    Mat.e=build(d.e-2)       #build
 	Index1.e=rep(0,ntri.e*Dsq.e)
 	Index2.e=Index1.e
 	S.e=Index1.e
 	place.e=1
#        browser()
 	for (k in 1:ntri.e){
 		LocK.e=locEng(V.e[T.e[k,1],],V.e[T.e[k,2],],V.e[T.e[k,3],],Mat.e,d.e,index.e)
 		result=which(LocK.e!=0,arr.in=TRUE)
 		i.e=result[,1];j.e=result[,2];s.e=LocK.e[cbind(i.e,j.e)]
  		L.e=length(i.e)
 		Index1.e[place.e:(place.e+L.e-1)]=(k-1)*D.e+i.e
 		Index2.e[place.e:(place.e+L.e-1)]=(k-1)*D.e+j.e
 		S.e[place.e:(place.e+L.e-1)]=s.e
 		place.e=place.e+L.e
 		}
 		K.e=matrix(0,ntri.e*D.e,ntri.e*D.e)
 		K.e[cbind(Index1.e[1:(place.e-1)],Index2.e[1:(place.e-1)])]=S.e[1:(place.e-1)]
 		return(K.e)
 	}


######################## locEng
########################
locEng=function(v1.lo,v2.lo,v3.lo,Mat.lo,d.lo,index.lo){
	D.lo=(d.lo+1)*(d.lo+2)/2
	Id.lo=diag(rep(1,D.lo))
	vx.lo=c(1,0)
	vy.lo=c(0,1)
    result=tcord(v1.lo,v2.lo,v3.lo,vx.lo)       #tcord
	lam1x.lo=result[,1];lam2x.lo=result[,2];lam3x.lo=result[,3]
	result=tcord(v1.lo,v2.lo,v3.lo,vy.lo)
	lam1y.lo=result[,1];lam2y.lo=result[,2];lam3y.lo=result[,3]
	Dx.lo=dirder(Id.lo,lam1x.lo,lam2x.lo,lam3x.lo)
	Dxx.lo=dirder(Dx.lo,lam1x.lo,lam2x.lo,lam3x.lo)
	Dxy.lo=dirder(Dx.lo,lam1y.lo,lam2y.lo,lam3y.lo)
	Dy.lo=dirder(Id.lo,lam1y.lo,lam2y.lo,lam3y.lo)
	Dyy.lo=dirder(Dy.lo,lam1y.lo,lam2y.lo,lam3y.lo)
#        browser()
	if(index.lo==1){
	K3.lo=abs(triarea(v1.lo,v2.lo,v3.lo))*(t(Dxx.lo)%*%Mat.lo%*%Dxx.lo+2*t(Dxy.lo)%*%Mat.lo%*%Dxy.lo+t(Dyy.lo)%*% Mat.lo%*%Dyy.lo)
	}
	if(index.lo==2){
	K3.lo=abs(triarea(v1.lo,v2.lo,v3.lo))*(t(Dxx.lo+Dyy.lo)%*%Mat.lo%*%(Dxx.lo+Dyy.lo))
	}
	return(K3.lo)
	}

######################## tcord
######################## return the barycentric coord of direction
tcord=function(v1.tc,v2.tc,v3.tc,v.tc){
	result1=bary(v1.tc,v2.tc,v3.tc,v.tc[1],v.tc[2])
	result2=bary(v1.tc,v2.tc,v3.tc,0,0)
	return(result1-result2)
	}

######################## build
######################## evaluate the matrix for inner product
build=function(d.bu){
	result=indices(d.bu)
	I.bu=result[,1];J.bu=result[,2];K.bu=result[,3]
	m.bu=(d.bu+1)*(d.bu+2)/2
	Mat.bu=matrix(0,m.bu,m.bu)
	for (j in 1:m.bu){
		for(k in 1:m.bu){
			Mat.bu[k,j]=choose(I.bu[j]+I.bu[k],I.bu[j])*choose(J.bu[j]+J.bu[k],J.bu[j])*choose(K.bu[j]+K.bu[k],K.bu[j])
			}
		}
	Mat.bu=Mat.bu/(choose(d.bu*2,d.bu)*choose(2*d.bu+2,2))
	return(Mat.bu)
	}


####################### dirder
#######################
dirder=function(Bc.di,lam1.di,lam2.di,lam3.di){
	m.di=dim(Bc.di)[1]
	d.di=degree(m.di)
	DerBc.di=d.di*de.cast.step(lam1.di,lam2.di,lam3.di,Bc.di)
	return(DerBc.di)
	}


################## refine
##################
refine=function(V0.re,T0.re,index.re,type.re){
	if(index.re=="total") index.re=1:dim(T0.re)[1]
	n.re=nrow(V0.re)
	m.re=nrow(T0.re)
	V1.re=V0.re
	T1.re=T0.re

	if(type.re=="four"){
	for (j in index.re){
		J.re=T0.re[j,]
		x1.re=V0.re[J.re[1],1];y1.re=V0.re[J.re[1],2]
		x2.re=V0.re[J.re[2],1];y2.re=V0.re[J.re[2],2]
		x3.re=V0.re[J.re[3],1];y3.re=V0.re[J.re[3],2]
		x12.re=(x1.re+x2.re)/2;y12.re=(y1.re+y2.re)/2
		x13.re=(x1.re+x3.re)/2;y13.re=(y1.re+y3.re)/2
		x23.re=(x2.re+x3.re)/2;y23.re=(y2.re+y3.re)/2

		place12=which((V1.re[,1]==x12.re)*(V1.re[,2]==y12.re)==1)
		place13=which((V1.re[,1]==x13.re)*(V1.re[,2]==y13.re)==1)
		place23=which((V1.re[,1]==x23.re)*(V1.re[,2]==y23.re)==1)
		if (length(place12)==0) {
			n.re=n.re+1
			place12=n.re
			V1.re=rbind(V1.re,c(x12.re,y12.re))
			}
		if (length(place13)==0){
			n.re=n.re+1
			place13=n.re
			V1.re=rbind(V1.re,c(x13.re,y13.re))
			}
		if (length(place23)==0){
			n.re=n.re+1
			place23=n.re
			V1.re=rbind(V1.re,c(x23.re,y23.re))
			}
		T1.re=rbind(T1.re,c(J.re[1],place12,place13))
		T1.re=rbind(T1.re,c(place12,place13,place23))
		T1.re=rbind(T1.re,c(J.re[2],place12,place23))
		T1.re=rbind(T1.re,c(J.re[3],place13,place23))
		}
		}
		if(type.re=="three"){
		for(j in index.re){
			J.re=T0.re[j,]
			x1.re=V0.re[J.re[1],1];y1.re=V0.re[J.re[1],2]
			x2.re=V0.re[J.re[2],1];y2.re=V0.re[J.re[2],2]
			x3.re=V0.re[J.re[3],1];y3.re=V0.re[J.re[3],2]
			x4.re=(x1.re+x2.re+x3.re)/3;y4.re=mean(y1.re+y2.re+y3.re)/3
			place=which((V1.re[,1]==x4.re)*(V1.re[,2]==y4.re)==1)
			if (length(place)==0){
			n.re=n.re+1
			place=n.re
			V1.re=rbind(V1.re,c(x4.re,y4.re))
			}
			T1.re=rbind(T1.re,c(J.re[1],J.re[2],place))
			T1.re=rbind(T1.re,c(J.re[2],J.re[3],place))
			T1.re=rbind(T1.re,c(J.re[3],J.re[1],place))
		}
		}
		#T1.re=T1.re[-1,]
		T1.re=T1.re[-index.re,]
		T1.re=mkcc(V1.re,T1.re)
		return(list(V1.re,T1.re))
	}

################# mkcc
#################
mkcc=function(V.mk,T.mk){
	m.mk=nrow(T.mk)
	for (j in 1:m.mk){
		if (triarea(V.mk[T.mk[j,1],],V.mk[T.mk[j,2],],V.mk[T.mk[j,3],])<0){
			T.mk[j,]=c(T.mk[j,1],T.mk[j,3],T.mk[j,2])
			}
		}
		return(T.mk)
	}


