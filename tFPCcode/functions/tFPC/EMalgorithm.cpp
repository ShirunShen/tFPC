#include<RcppArmadillo.h>
#include<Rcpp.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
arma::mat Tcreate(const arma::vec& K0_Tc,
                  const int& J_Tc,
                  const int& p_Tc){
    //use for defining T matrix of state equation
    arma::mat temp = kron(K0_Tc.t(), eye(J_Tc,J_Tc));
    temp = join_cols(temp, eye(p_Tc*J_Tc,p_Tc*J_Tc));
    temp = join_rows(temp, zeros((p_Tc+1)*J_Tc,J_Tc));
    return(temp);
}

//[[Rcpp::export]]
arma::mat sub_m(const arma::mat& A_m, const arma::vec& n_m, const int& i){
    
    arma::mat B_m;
    if(i == 1) B_m = A_m.rows(0, n_m(0)-1);
    else B_m = A_m.rows(sum(n_m.subvec(0,i-2)), sum(n_m.subvec(0,i-1))-1);
    
    return B_m;
} //output the sub_rows of matrix

//select subvector of z, etc.

//[[Rcpp::export]]
arma::vec sub_v(const arma::vec& A_v, const arma::vec& n_v, const int& i){
    
    arma::vec B_v;
    if(i == 1) B_v = A_v.rows(0, n_v(0)-1);
    else B_v = A_v.rows(sum(n_v.subvec(0,i-2)), sum(n_v.subvec(0,i-1))-1);
    
    return B_v;
}

//[[Rcpp::export]]
List kalmanfilter(const arma::vec& z_kal,
                  const arma::mat& B_kal,
                  const arma::vec& theta_b_kal,
                  const arma::vec& theta_c_kal,
                  const arma::mat& c_kal,
                  const arma::vec& b1_kal,
                  const arma::mat& Q1_kal,
                  const arma::mat& HJ0_kal,
                  const arma::mat& Theta0_kal,
                  const arma::mat& T_kal,
                  const arma::vec& ni_kal,
                  const double sigma20_kal){
    int n_kal = ni_kal.n_elem;
    int J_kal = HJ0_kal.n_rows;
    int p_kal = T_kal.n_rows/J_kal - 1;
    int ntmax = ni_kal.max();
    int count = 0;
    arma::vec bt = b1_kal;
    arma::mat Qt = Q1_kal;
    arma::vec b_kal  = bt;
    arma::mat Q_kal = Qt;
    arma::mat invF_kal = zeros(n_kal * ntmax, ntmax);
    arma::vec v_kal = zeros(sum(ni_kal));
    arma::mat Theta_bc_kal = theta_b_kal * theta_c_kal.t();
    for(int t=1; t<= n_kal; t++){
        
        arma::vec zt = sub_v(z_kal, ni_kal, t);
        arma::mat ct = c_kal.row(t-1).t(); //.row returns matrix
        arma::mat Bt = sub_m(B_kal, ni_kal, t);
        int       nt = ni_kal(t-1);
        arma::mat temp = zeros(nt, p_kal * J_kal);
        arma::mat BtTh = Bt * Theta0_kal;
        arma::mat BtThS = join_rows(temp,BtTh);
        
        arma::vec vt = zt - Bt * Theta_bc_kal * ct - BtThS * bt;
        arma::mat Ft = BtTh * Qt.submat(p_kal*J_kal,p_kal*J_kal,
                (p_kal+1)*J_kal-1,(p_kal+1)*J_kal-1) * BtTh.t() + sigma20_kal
                * eye(nt,nt);
        arma::mat invFt = inv(Ft);
        arma::vec btt = bt + Qt * BtThS.t() * invFt * vt;
        arma::mat Qtt = Qt - Qt * BtThS.t() * invFt * BtThS * Qt;
        
        bt = T_kal * btt;
        Qt = T_kal * Qtt * T_kal.t(); //with following for loop
        for(int j=1; j<=J_kal; j++){
            Qt(j-1,j-1) = Qt(j-1,j-1) + HJ0_kal(j-1,j-1);
        }
        
        b_kal = join_cols(b_kal, bt);
        Q_kal = join_cols(Q_kal, Qt);
        invF_kal.submat((t-1)*ntmax, 0, (t-1)*ntmax+nt-1, nt-1) = invFt;
        v_kal.subvec(count,count+nt-1) = vt;
        count = count + nt;
    }
    
    List res = List::create(Named("b") = b_kal,
                            Named("Q") = Q_kal,
                            Named("invF") = invF_kal,
                            Named("v") = v_kal);
    
    return res;
}


//[[Rcpp::export]]
List statesmooth(const arma::mat& B_sta,
                 const arma::mat& Theta0_sta,
                 const arma::mat& invF_sta,
                 const arma::mat& Q_sta,
                 const arma::vec& b_sta,
                 const arma::vec& v_sta,
                 const arma::vec& ni_sta,
                 const arma::mat& T_sta){
    int n_sta = ni_sta.n_elem;
    int J = Theta0_sta.n_cols;
    int p = T_sta.n_rows/J - 1;
    int ntmax = ni_sta.max();
    arma::vec rt = zeros(T_sta.n_cols);
    arma::mat Nt = zeros(T_sta.n_cols,T_sta.n_cols);
    arma::vec bhat  = zeros(n_sta*(p+1)*J);
    arma::mat Vhat  = zeros(n_sta*J,(p+1)*J);
    arma::vec alphahat = zeros(n_sta*J);
    arma::mat Sigmahat = zeros(n_sta*J, J);
    
    for(int t=n_sta; t>=1; t--){
        int       nt   = ni_sta(t-1);
        arma::mat Bt   = sub_m(B_sta,ni_sta,t);
        arma::mat temp = zeros(nt, p*J);
        arma::mat BtTh = Bt * Theta0_sta;
        arma::mat BtThS = join_rows(temp,BtTh);
        arma::mat invFt = invF_sta.submat((t-1)*ntmax,0,(t-1)*ntmax+nt-1, nt-1);
        arma::mat Qt    = Q_sta.rows((t-1)*(p+1)*J, t*(p+1)*J-1);
        arma::vec bt    = b_sta.subvec((t-1)*(p+1)*J, t*(p+1)*J-1);
        arma::vec vt    = sub_v(v_sta,ni_sta,t);
        arma::mat Lt    = T_sta - T_sta * Qt * BtThS.t() * invFt * BtThS;
        
        rt = BtThS.t() * invFt * vt + Lt.t() * rt;
        Nt = BtThS.t() * invFt * BtThS + Lt.t() * Nt * Lt;
        arma::vec bthat = bt + Qt * rt;
        arma::mat Vthat = Qt - Qt * Nt * Qt;
        bhat.subvec((t-1)*(p+1)*J, t*(p+1)*J-1) = bthat;
        Vhat.submat((t-1)*J,0,t*J-1,(p+1)*J-1) = Vthat.submat(p*J,0,(p+1)*J-1,(p+1)*J-1);
        alphahat.subvec((t-1)*J,t*J-1) = bthat.subvec(p*J, (p+1)*J-1);
        Sigmahat.submat((t-1)*J,0,t*J-1,J-1) = Vthat.submat(p*J,p*J,(p+1)*J-1,(p+1)*J-1);
    } 
    
    List res = List::create(Named("alphahat") = alphahat,
                            Named("Sigmahat") = Sigmahat,
                            Named("Vhat")     = Vhat);
    
    return res;
}





//[[Rcpp::export]]
double obj_func(const arma::vec& x,
                const arma::vec& m,
                const arma::mat& A){
    arma::mat temp = (x-m).t()*A*(x-m);
    return temp(0,0);
}

//[[Rcpp::export]]
arma::vec optim_sphere(const arma::vec& thetac0_o,
                       const arma::vec& m_o,
                       const arma::mat& A_o){
    arma::vec theta_c_o = thetac0_o; //initial value and the container of results
    
    double alpha_bar = 1;
    double beta = 0.1;
    double sigma = 0.1;
    
    for(int iter=0; iter<=200; iter++){
        arma::vec theta_ck = theta_c_o;
        arma::vec eta_k = - 2 * (A_o * (theta_ck - m_o) 
                                     - theta_ck * theta_ck.t() * A_o * (theta_ck - m_o));
        
        double criteria = -1;
        int           n = 0;
        int    stopping = 1;
        arma::vec Rxk = zeros(size(theta_ck));
        while(criteria < -1e-10 && stopping < 100){
            double coef = alpha_bar * pow(beta,n);
            Rxk = theta_ck + alpha_bar * pow(beta,n) * eta_k;
            
            double fRxk = obj_func(Rxk/norm(Rxk),m_o,A_o);
            double fxk = obj_func(theta_ck,m_o,A_o);
            criteria = fxk - fRxk - sigma * coef * sum(eta_k % eta_k);
            n++;    
            stopping++;
            //cout << n << endl;
            //cout << criteria << endl;
        } //loop while
        theta_c_o = Rxk/norm(Rxk);
        //cout << iter << endl;
        //cout << "diff" << norm(theta_c_o - theta_ck) << endl;
        
        if(norm(theta_c_o - theta_ck) < 0.01){ 
            //cout << "the number of iteration for thetac" << iter << endl;
            break;
        }
        
    } //loop for
    
    return theta_c_o;
} //Armijo line search for the Rayleigh quotient on sphere


//[[Rcpp::export]]
arma::vec theta_b_update(const arma::vec& z_t,
                         const arma::mat& B_t,
                         const arma::mat& c_t,
                         const arma::vec& theta_c_t,
                         const arma::vec& thetab0_t,
                         const arma::mat& Theta_t,
                         const arma::vec& alphahat_t,
                         const arma::mat& Ps_t,
                         const arma::vec& ni_t,
                         const double& sigma2_t,
                         const double& lambmu1_t){
    int n_t = ni_t.n_elem;
    int J = Theta_t.n_cols;
    int nb = B_t.n_cols;
    //int nc = c_t.n_cols;
    arma::vec temp = zeros(nb);
    arma::mat temp2 = zeros(nb, nb);
    arma::vec theta_b_t = zeros(nb);
    
    
    for(int t = 1; t <= n_t; t++){
        arma::vec zt = sub_v(z_t, ni_t, t);
        arma::mat Bt = sub_m(B_t, ni_t, t);
        arma::vec ct = c_t.row(t-1).t();
        arma::vec alphahatt = alphahat_t.subvec((t-1)*J, t*J-1);
        arma::mat tem = theta_c_t.t() * ct;
        double    temp3 = tem(0,0);
        temp = temp + temp3 * Bt.t() * (zt - Bt * Theta_t * alphahatt);
        temp2 = temp2 + (temp3 * temp3) * Bt.t() * Bt;
    }
    
    temp2 = temp2 + sigma2_t * lambmu1_t * Ps_t;
    
    //test
    //cout << "range of component 1" << max(temp2) << min(temp2) << endl;
    //cout << "range of penalty matrix" << max(Ps_t) << min(Ps_t) << endl;
    arma::mat A_t = temp2;
    arma::vec m_t = inv(A_t) * temp;
    
    theta_b_t = optim_sphere(thetab0_t,m_t,A_t); //inv(temp2) * temp;
    //theta_b_t = solve(temp2,temp);
    return theta_b_t;
}





//[[Rcpp::export]]
arma::vec theta_c_update(const arma::vec& z_t,
                         const arma::mat& B_t,
                         const arma::mat& c_t,
                         const arma::vec& theta_b_t,
                         const arma::mat& Theta_t,
                         const arma::vec& alphahat_t,
                         const arma::mat& Pt_t,
                         const arma::vec& ni_t,
                         const double& sigma2_t,
                         const double& lambmu2_t){
    int n_t = ni_t.n_elem;
    int J = Theta_t.n_cols;
    int nc = c_t.n_cols;
    arma::vec temp = zeros(nc);
    arma::mat temp2 = zeros(nc,nc);
    arma::vec theta_c_t = zeros(nc);
    arma::vec m_t = zeros(nc);
    arma::mat A_t = zeros(nc,nc);
    for(int t=1; t<=n_t; t++){
        arma::vec zt = sub_v(z_t, ni_t, t);
        arma::mat Bt = sub_m(B_t, ni_t, t);
        arma::vec ct = c_t.row(t-1).t();
        arma::vec alphahatt = alphahat_t.subvec((t-1)*J, t*J-1);
        
        arma::mat temp3 = Bt * theta_b_t * ct.t();
        temp = temp + temp3.t() * (zt - Bt * Theta_t * alphahatt);
        temp2 = temp2 + temp3.t() * temp3;
    } //loop for
    A_t = temp2 + sigma2_t * lambmu2_t * Pt_t;
    m_t = inv(A_t) * temp;
    // the following is optimization on sphere
    theta_c_t = m_t; //optim_sphere(thetac0_t,m_t,A_t);
    return theta_c_t;
}



//[[Rcpp::export]]
double sigma2_update(const arma::vec& z_s,
                     const arma::vec& thetab_s,
                     const arma::vec& thetac_s,
                     const arma::mat& B_s,
                     const arma::mat& c_t_s,
                     const arma::mat& Theta_s,
                     const arma::vec& alphahat_s,
                     const arma::mat& Sigmahat_s,
                     const arma::vec& ni_s){
    
    double sigma2hat = 0;
    int    J   = Theta_s.n_cols;
    int    n_s = ni_s.n_elem;
    arma::mat Theta_bc_s = thetab_s * thetac_s.t();
    
    for(int t=1;t<=n_s;t++){
        arma::vec zt = sub_v(z_s, ni_s, t);
        arma::mat Bt = sub_m(B_s, ni_s, t);
        arma::vec ct = c_t_s.row(t-1).t();
        arma::vec alphat = alphahat_s.subvec((t-1)*J, t*J-1);
        arma::mat Sigmat = Sigmahat_s.submat((t-1)*J, 0, t*J-1, J-1);
        
        arma::vec temp  = zt - Bt * Theta_bc_s * ct - Bt * Theta_s * alphat;
        temp = temp.t() * temp + trace(Bt * Theta_s * Sigmat * Theta_s.t() * Bt.t());
        sigma2hat = sigma2hat + temp(0); //change 1*1 vector into double:(0)
    }
    sigma2hat = sigma2hat/sum(ni_s);
    return sigma2hat;
}


//[[Rcpp::export]]
arma::mat Theta_update(const arma::vec& z_th,
                       const arma::mat& c_th,
                       const arma::mat& B_th,
                       const arma::vec& theta_b_th,
                       const arma::vec& theta_c_th,
                       const arma::mat& Theta_th,
                       const arma::vec& alphahat_th,
                       const arma::mat& Sigmahat_th,
                       const arma::mat& P_th,
                       const arma::vec& ni_th,
                       const double sigma2_th,
                       const double lambpc_th){
    int n_th = ni_th.n_elem;
    int J_th = Theta_th.n_cols;
    int nb_th = B_th.n_cols;
    arma::mat Theta_bc_th = theta_b_th * theta_c_th.t();
    arma::mat Thetahat_th = zeros(Theta_th.n_rows, Theta_th.n_cols);
    for(int j=1; j<=J_th; j++){
        arma::mat temp = zeros(nb_th,nb_th);
        arma::vec temp2 = zeros(nb_th);
        for(int t=1; t<=n_th; t++){
            arma::mat Bt = sub_m(B_th, ni_th, t);
            arma::vec zt = sub_v(z_th, ni_th, t);
            arma::vec ct = c_th.row(t-1).t();
            arma::vec alphat = alphahat_th.subvec((t-1)*J_th, t*J_th-1);
            arma::mat Sigmat = Sigmahat_th.submat((t-1)*J_th, 0, t*J_th-1, J_th-1);
            
            temp = temp + (alphat(j-1)*alphat(j-1) + Sigmat(j-1,j-1)) * Bt.t() * Bt;
            arma::vec temp3 = zeros(ni_th(t-1));
            for(int k=1; (k!=j) && k <= J_th; k++){
                temp3 = temp3 + (alphat(k-1)*alphat(j-1) + Sigmat(k-1,j-1)) * Bt * Theta_th.col(k-1);
            }  //k iteration
            
            temp2 = temp2 + Bt.t() * ((zt - Bt * Theta_bc_th * ct) * alphat(j-1) - temp3);
        }  //t iteration
        Thetahat_th.col(j-1) = inv(temp + sigma2_th * lambpc_th * P_th) * temp2;
        //Thetahat_th.col(j-1) = solve(temp+sigma2_th*lambpc_th*P_th, temp2);
    }  //j iteration
    return Thetahat_th;
}







double D_elem(const int& i,
              const int& k,
              const int& j,
              const int& n_D,
              const arma::vec& alphahat_D,
              const arma::mat& Vhat_D){
    double elem = 0.0;
    int    J    = alphahat_D.n_elem/n_D;
    int    p    = Vhat_D.n_cols/J - 1;

    if(k-i>=0){
        for(int l=i; l<=n_D + 1 - k; l++){
            elem = elem + alphahat_D((l-1)*J+j-1) * alphahat_D((k+l-i-1)*J+j-1) + Vhat_D((l-1)*J+j-1,(p-k+i)*J+j-1);
        }
    }else{
        for(int l=i; l<=n_D + 1 - k; l++){
            elem = elem + alphahat_D((l-1)*J+j-1) * alphahat_D((k+l-i-1)*J+j-1) + Vhat_D((l+k-i-1)*J+j-1,(p-i+k)*J+j-1); //need to revise here.
        }
    }
    
    return elem;
}


arma::mat Dhat_p(const int& j,
                 const int& n_D,
                 const int& p,
                 const arma::vec& alphahat_D,
                 const arma::mat& Vhat_D){
    arma::mat D = zeros(p,p);
    for(int m=2; m<=p+1; m++){
        for(int n=2; n<=p+1; n++){
            D(m-2,n-2) = D_elem(m,n,j,n_D, alphahat_D,Vhat_D);
        }
    }
    return D;
}

arma::vec dhat_j(const int& j,
                 const int& n_D,
                 const int& p,
                 const arma::vec& alphahat_d,
                 const arma::mat& Vhat_d){
    arma::vec elem = zeros(p);
    for(int q = 2; q<=p+1; q++){
        elem(q-2) = D_elem(1, q, j, n_D, alphahat_d,Vhat_d);
    }
    return elem;
}




//for fliping the columns of matrix

arma::mat flip_m(const arma::mat& A_fm){
    int ncol = A_fm.n_cols;
    arma::mat B_fm = zeros(A_fm.n_rows, A_fm.n_cols);
    for (int i = 1; i <= ncol; i++){
        B_fm.col(ncol-i) =  A_fm.col(i-1);
    }
    return B_fm;
}

arma::vec flip_v(const arma::vec& A_fv){
    int num = A_fv.n_elem;
    arma::vec B_fv = zeros(A_fv.n_elem);
    for (int i = 1; i <= num; i++){
        B_fv(num-i) = A_fv(i-1);
    }
    return B_fv;
}

//[[Rcpp::export]]
arma::vec theta_b_init(const arma::vec& z_t,
                       const arma::mat& B_t,
                       const arma::mat& c_t,
                       const arma::vec& thetab0_t,
                       const arma::vec& theta_c_t,
                       const arma::vec& ni_t,
                       const arma::mat& Ps_t,
                       const double lambda){
    int n_t = ni_t.n_elem;
    int nb  = B_t.n_cols;
    
    arma::vec temp = zeros(nb);
    arma::mat temp2 = zeros(nb,nb);
    arma::vec theta_b_t = zeros(nb);
    
    for(int t=1; t<=n_t; t++){
        arma::vec zt = sub_v(z_t, ni_t, t);
        arma::mat Bt = sub_m(B_t, ni_t, t);
        arma::vec ct = c_t.row(t-1).t();
        arma::mat tem = theta_c_t.t() * ct;
        double temp3 = tem(0,0);
        temp = temp + temp3 * Bt.t() * zt;
        temp2 = temp2 + (temp3 * temp3) * Bt.t() * Bt;
    }
    temp2 = temp2 + lambda * Ps_t;
    arma::mat A_t = temp2;
    arma::vec m_t = inv(A_t) * temp;
    
    //theta_b_t = m_t;
    theta_b_t = optim_sphere(thetab0_t,m_t,A_t);
    
    return theta_b_t;
}

//[[Rcpp::export]]
arma::vec theta_c_init(const arma::vec& z_t,
                       const arma::mat& B_t,
                       const arma::mat& c_t,
                       const arma::vec& theta_b_t,
                       const arma::vec& thetac0_t,
                       const arma::vec& ni_t,
                       const arma::mat& Pt_t,
                       const double lambda){
    int n_t = ni_t.n_elem;
    int nc = c_t.n_cols;
    arma::vec temp = zeros(nc);
    arma::mat temp2 = zeros(nc,nc);
    arma::vec theta_c_t = zeros(nc);
    
    for(int t=1; t<=n_t; t++){
        arma::vec zt = sub_v(z_t,ni_t,t);
        arma::mat Bt = sub_m(B_t,ni_t,t);
        arma::vec ct = c_t.row(t-1).t();
        
        arma::mat temp3 = Bt * theta_b_t * ct.t();
        temp = temp + temp3.t() * zt;
        temp2 = temp2 + temp3.t() * temp3;
    }
    temp2  = temp2 + lambda * Pt_t;
    arma::mat A_t = temp2;
    arma::vec m_t = inv(A_t) * temp;
    //theta_c_t = optim_sphere(thetac0_t,m_t,A_t);
    theta_c_t = m_t;
    
    return theta_c_t;
}



                       
                       
                       
                       
                       
                       
//[[Rcpp::export]]
List EMinit(const arma::vec& z_em,
            const arma::mat& B_em,
            const arma::mat& c_em,
            const arma::vec& ni_em,
            const arma::vec& thetab0_em,
            const arma::vec& thetac0_em,
            const arma::mat& Ps_em,
            const arma::mat& Pt_em,
            const double lambdab_em,
            const double lambdac_em){
    int n_em = ni_em.n_elem;
    int nbasis_em = B_em.n_cols;
    arma::vec thetab_em = thetab0_em;
    arma::vec thetac_em = thetac0_em;
    int Maxiter = 400;
    int iter = 1;
    
    double criteria = 1.0;
    
    while(iter <= Maxiter && criteria>1e-3){
        arma::vec thetab_old = thetab_em;
        arma::vec thetac_old = thetac_em;
        
        thetac_em = theta_c_init(z_em,B_em,c_em,thetab_em,thetac_em,ni_em,Pt_em,lambdac_em);
        thetab_em = theta_b_init(z_em,B_em,c_em,thetab_em,thetac_em,ni_em,Ps_em,lambdab_em);
        
        arma::vec criteria0 = zeros(2);
        criteria0(0) = norm(thetab_em - thetab_old)/norm(thetab_old);
        criteria0(1) = norm(thetac_em - thetac_old)/norm(thetac_old);
        criteria = max(criteria0);
        iter++;
    }
    //cout << iter << endl;
    return List::create(Named("thetab") = thetab_em,
                        Named("thetac") = thetac_em);
}



//[[Rcpp::export]]
List EMalgorithm(const arma::vec& z_em,
                 const arma::mat& B_em,
                 const arma::mat& c_em,
                 const arma::mat& Ps_em,
                 const arma::mat& Pt_em,
                 const arma::vec& ni_em,
                 const arma::mat& HJ0_em,
                 const arma::vec& thetab0_em,
                 const arma::vec& thetac0_em,
                 const arma::mat& Theta0_em,
                 const arma::vec& K0_em,
                 const double& sigma20_em,
                 const double& lambmus_em,
                 const double& lambmut_em,
                 const double& lambpc_em){
    
    int n_em = ni_em.n_elem;
    int J_em = HJ0_em.n_rows;
    int p_em = K0_em.n_elem;
    int nbasis_em = B_em.n_cols;
    
    double    sigma2_em = sigma20_em;
    arma::mat     HJ_em = HJ0_em;
    arma::vec thetab_em = thetab0_em;
    arma::vec thetac_em = thetac0_em;
    arma::mat  Theta_em = Theta0_em;
    arma::vec      K_em = K0_em;
    
    
    int Maxiter = 500;
    int iter    = 1;

    arma::vec seq_sigma = zeros(Maxiter);
    arma::vec seq_K     = zeros(p_em * Maxiter);
    arma::mat seq_HJ    = zeros(J_em * Maxiter, J_em);
    
    
    arma::vec b1_em = zeros((p_em+1)*J_em);
    arma::mat Q1_em = 0.1 * eye((p_em+1)*J_em,(p_em+1)*J_em);
    
    arma::vec alphahat_em = zeros(J_em * n_em);
    arma::mat Sigmahat_em = zeros(J_em * n_em, J_em);
    arma::mat Vhat_em     = zeros(J_em * n_em, J_em * (p_em+1));
    double criteria = 1.0;
    
    while(iter <= Maxiter && criteria > 1e-2){ //Sep 12, change 1e-3 to 1e-2
   
        double sigma02_em = sigma2_em;
        
        arma::mat T_em = Tcreate(K_em, J_em, p_em);
        //cout << T_em << endl;

        
        //kalman filter
        List filter = kalmanfilter(z_em,B_em,thetab_em,thetac_em,c_em,b1_em,
                      Q1_em,HJ_em,Theta_em,T_em,ni_em,sigma2_em);

        arma::vec b_em = filter["b"];
        //cout << max(b_em) << min(b_em) << endl;
        arma::mat Q_em = filter["Q"];
        //cout << max(Q_em) << min(Q_em) << endl;
        arma::mat invF_em = filter["invF"];
        //cout << invF_em << endl;
        arma::vec v_em = filter["v"];
        //cout << max(v_em) << min(v_em) << endl;
        
        
        //kalman smoother
        List smooth = statesmooth(B_em,Theta_em,invF_em,
                                  Q_em,b_em,v_em,ni_em,T_em);
        arma::vec alphahat0_em = smooth["alphahat"];
        alphahat_em = alphahat0_em;
        arma::mat Sigmahat0_em = smooth["Sigmahat"];
        Sigmahat_em = Sigmahat0_em;
        arma::mat Vhat0_em = smooth["Vhat"];
        Vhat_em     = Vhat0_em;   //full covariance matrix with dim n*J, (p+1)*J
        //cout << "Vhat_em" << Vhat_em << endl;
        
        
        
        //to make sure the expectations of alpha equal to zero.
        arma::vec temp_ave = zeros(J_em);
        for(int j=1; j<=J_em; j++){
            for(int l=1; l<=n_em; l++){
               temp_ave(j-1) = temp_ave(j-1) + alphahat_em(J_em*(l-1)+(j-1));
            }
        }
        temp_ave = temp_ave/n_em;
        
        for(int j=1; j<=J_em; j++){
            for(int l=1; l<=n_em; l++){
        //        alphahat_em(J_em*(l-1)+(j-1)) = alphahat_em(J_em*(l-1)+(j-1)) -  temp_ave(j-1);
            }
        }
        //cout << min(temp_ave) << max(temp_ave) << endl;

        
        //cout << mean(alphahat_em) << endl;
        //cout << max(Sigmahat_em) << endl;
        int iter_Mstep = 1;
        double criteria_Mstep = 1.0;
        while(iter_Mstep <= 100 && criteria_Mstep > 1e-2){ //Apr 30, change 1e-1 to 1e-2
            
            double sigma02Mstep_em = sigma2_em;
            arma::vec thetabMstep_em = thetab_em;
            arma::vec thetacMstep_em = thetac_em;
            arma::mat ThetaMstep_em  = Theta_em;
            
            
            //update thetac
            thetac_em = theta_c_update(z_em,B_em,c_em,thetab_em,Theta_em,
                                       alphahat_em,Pt_em,ni_em,sigma2_em,lambmut_em); //normalization
            //cout << min(thetac_em) << " " << max(thetac_em) << endl;
            
            
            //update thetab
            thetab_em = theta_b_update(z_em,B_em,c_em,thetac_em,thetab_em, Theta_em,alphahat_em,Ps_em,ni_em,sigma2_em,lambmus_em);
            //cout << thetab_em << endl;
            
            
            //update sigma2
            sigma2_em = sigma2_update(z_em,thetab_em,thetac_em,B_em,c_em,Theta_em,
                                      alphahat_em,Sigmahat_em,ni_em);
            //cout << sigma2_em << endl;
            
            
            //update Theta without normalization
            Theta_em = Theta_update(z_em,c_em,B_em,thetab_em,thetac_em,Theta_em,
                alphahat_em,Sigmahat_em,Ps_em,ni_em,sigma2_em,lambpc_em);
            //cout << Theta_em << endl;
            
            
            
            //update HJ without normalization
            for(int j=1; j<=J_em; j++){
                arma::mat Dhat_em = zeros(p_em+1, p_em+1);
                
                for(int m=1; m<=p_em+1; m++){
                    for(int n=1; n<=p_em+1; n++){
                        Dhat_em(m-1,n-1) = D_elem(m,n,j,n_em,alphahat_em, Vhat_em);
                        if(m==1 && n!=1){
                            Dhat_em(m-1,n-1) = - Dhat_em(m-1,n-1);
                        }
                        if(m!=1 && n==1){
                            Dhat_em(m-1,n-1) = - Dhat_em(m-1,n-1);
                        }
                    }
                }
                //cout << Dhat_em << endl;
                //for(int m=1; m<=p_em+1; m++){
                //    for(int n=1; n<=p_em+1; n++){
                //        Dhat_em(m-1,n-1) = D_elem(m,n,j,n_em,alphahat_em, Vhat_em);
                //    }
                //}
                arma::vec K_tilde = ones(p_em+1);
                for(int i=1; i<=p_em; i++) K_tilde(i) = K_em(i-1);
                arma::mat S_j = K_tilde.t() * Dhat_em * K_tilde;
                HJ_em(j-1,j-1) = S_j(0)/n_em;
            }
            //cout << HJ_em << endl;
            
            
            
            //update K
            arma::mat Dphat_em = zeros(p_em, p_em);
            arma::vec dhat_em = zeros(p_em);
            
            for(int j=1; j<= J_em; j++){
                //Here is the test result, keep watching!
                //Dphat_em = Dphat_em + Dhat_p(j, n_em, p_em, alphahat_em);
                //dhat_em = dhat_em + dhat_j(j, n_em, p_em, alphahat_em);
                Dphat_em = Dphat_em + Dhat_p(j, n_em, p_em, alphahat_em, Vhat_em)/HJ_em(j-1,j-1);
                
                dhat_em = dhat_em + dhat_j(j, n_em, p_em, alphahat_em, Vhat_em)/HJ_em(j-1,j-1);
                
            }
            K_em = inv(Dphat_em) * dhat_em;
            //K_em = solve(Dphat_em, dhat_em);
            
            //cout << K_em << endl;
            
            
            //update HJ and Theta with normalization
            arma::vec eigvalue = zeros(nbasis_em);
            arma::mat eigvec = zeros(nbasis_em,nbasis_em);
            
            //test Mar 7, 2020
            arma::mat Thetatilde_em = Theta_em; //for updating alphahat_em;
            //end test
            
            eig_sym(eigvalue, eigvec, Theta_em * HJ_em * Theta_em.t());
            HJ_em = diagmat(flip_v(eigvalue.subvec(nbasis_em-J_em, nbasis_em - 1)));
            Theta_em = flip_m(eigvec.cols(nbasis_em - J_em, nbasis_em - 1));
            
            
            //test Mar 7, 2020, to transform 
            arma::mat trans_alpha_em = Theta_em.t() * Thetatilde_em;
            
            
            for(int t=1; t<= n_em; t++){
                alphahat_em.subvec((t-1)*J_em,t*J_em-1) = trans_alpha_em * alphahat_em.subvec((t-1)*J_em,t*J_em-1);
            }
            //end test
            
            
            //cout << Theta_em << endl;
            //cout << HJ_em << endl;
            
            arma::vec criteria0_Mstep = zeros(2);
            //criteria0_Mstep(0) = abs(sigma2_em - sigma02Mstep_em);
            criteria0_Mstep(0) = norm(thetabMstep_em - thetab_em)/norm(thetabMstep_em);
            criteria0_Mstep(1) = norm(thetacMstep_em - thetac_em)/norm(thetacMstep_em);
            //criteria0_Mstep(3) = norm(ThetaMstep_em + Theta_em)/norm(ThetaMstep_em);
            criteria_Mstep = max(criteria0_Mstep);
            
            //criteria_Mstep = 1;
            //criteria_Mstep = criteria0_Mstep(3);
            //cout << criteria_Mstep << endl;
            
            iter_Mstep++;
        }
        
        //test
        //cout << "sigma " << sigma2_em << endl;
        
        
        seq_sigma(iter-1) = sigma2_em;
        seq_K.subvec((iter-1)*p_em, iter*p_em -1) = K_em;
        seq_HJ.submat((iter-1)*J_em, 0, size(HJ_em)) = HJ_em;
        
        if(iter >=2){
            arma::vec criteria0=zeros(2);
            criteria0(0) = abs(sigma2_em - sigma02_em)/(abs(sigma02_em)+0.01);
            criteria0(1) = norm(K_em - seq_K.subvec((iter-2)*p_em, (iter-1)*p_em-1))/(norm(seq_K.subvec((iter-2)*p_em, (iter-1)*p_em-1))+0.01);
            criteria = max(criteria0);
        }   // for stopping criteria
        
        //cout << criteria << endl;
        //criteria = abs(sigma2_em - sigma02_em)/(abs(sigma02_em)+0.00001);
        
        //March 24, 2020, test for iteration convergence
        //criteria = 1.0;
        
        
        iter++;
    } //loop while
    
    
    return List::create(Named("seq_sigma")= seq_sigma,
                        Named("seq_K")    = seq_K,
                        Named("seq_HJ")   = seq_HJ,
                        Named("thetab")   = thetab_em,
                        Named("thetac")   = thetac_em,
                        Named("Theta")    = Theta_em,
                        Named("sigma2")   = sigma2_em,
                        Named("HJ")       = HJ_em,
                        Named("K")        = K_em,
                        Named("alphahat") = alphahat_em,
                        Named("Sigmahat") = Sigmahat_em,
                        Named("iter")     = iter-1);
}
