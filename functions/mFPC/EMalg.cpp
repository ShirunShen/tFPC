#include<RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::mat sub_m(const arma::mat& A_m, const arma::vec& n_m, const int& i);
arma::vec sub_v(const arma::vec& A_v, const arma::vec& n_v, const int& i);
arma::vec alpha_update(const arma::mat& D_al, const arma::mat& Theta_al, const arma::vec thetamu_al, const double sigma2_al, const arma::mat Bi_al, const arma::vec zi_al);
arma::mat Sigma_update(const arma::mat& D_si, const arma::mat& Theta_si, const double sigma2_si, const arma::mat& Bi_si);
arma::mat s2_update(const arma::mat& Theta_s, const arma::vec& thetamu_s, const arma::vec& zi_s, const arma::mat& Bi_s, const arma::vec& alphai_s, const arma::mat& Sigmai_s);
arma::mat flip_m(const arma::mat& A_fm);
arma::vec flip_v(const arma::vec& A_fv);


//[[Rcpp::export]]
Rcpp::List EM_algorithm(const arma::vec& z_em,
                        const arma::mat& B_em,
                        const arma::mat& P_em,
                        const arma::vec& ni_em,
                        const double& lambmu_em,
                        const double& lambpc_em,
                        const double& sigma20_em,
                        const arma::mat& D0_em,
                        const arma::vec& thetamu0_em,
                        const arma::mat& Theta0_em){
    
    int    n_em   = ni_em.n_elem;    
    int    J_em   = D0_em.n_rows;    
    int nbasis_em = B_em.n_cols;    
    
    double    sigma2_em  = sigma20_em;
    arma::mat D_em       = D0_em;
    arma::vec thetamu_em = thetamu0_em;
    arma::mat Theta_em   = Theta0_em;
    
    arma::mat alpha_em = zeros(J_em, n_em);  //store alpha
    arma::mat Sigma_em = zeros(J_em * n_em, J_em);
    arma::mat M_em     = zeros(J_em * n_em, J_em);
    
    double sigma2hat_em = 0.0;
    arma::mat Dhat_em  = zeros(J_em, J_em);
    arma::vec thetamuhat_em = zeros(nbasis_em);
    arma::mat Thetahat_em   = zeros(nbasis_em,J_em); //store Thetahat
    
    arma::vec eigvalue = zeros(nbasis_em);
    arma::mat eigvec = zeros(nbasis_em,nbasis_em);
    
    int iter = 1;
    arma::mat alpha0_em = ones(size(alpha_em));
    
    //major part of EM algorithm: iterations
    while(iter <= 120 && norm(alpha0_em-alpha_em)/(norm(alpha0_em)+0.01)>0.01){
    
        iter++;
        alpha0_em = alpha_em;
        
        
        double    temp  = 0.0;
        arma::mat Mi_em = zeros(J_em, J_em);
        arma::mat temp2 = zeros(J_em, J_em);
        arma::mat temp3 = zeros(nbasis_em, nbasis_em);
        arma::vec temp4 = zeros(nbasis_em);
        
        //test, June 7, 2021
        //arma::vec temp_ave = zeros(J_em);
        for(int i=1; i<=n_em; i++){
            arma::mat Bi_em     = sub_m(B_em, ni_em, i);
            arma::vec zi_em     = sub_v(z_em, ni_em, i);
            arma::vec alphai_em = alpha_update(D_em, Theta_em, thetamu_em, sigma2_em, Bi_em, zi_em);
            arma::mat Sigmai_em = Sigma_update(D_em, Theta_em, sigma2_em, Bi_em);
            
            //temp_ave = temp_ave + alphai_em;
            
            alpha_em.col(i-1)   = alphai_em;
            Sigma_em.rows(J_em*(i-1), J_em*i - 1) = Sigmai_em;
        }
        //temp_ave = temp_ave/n_em;
        
        //force alpha mean equals to zero.
        //for(int i=1; i<=n_em; i++){
            //alpha_em.col(i-1) = alpha_em.col(i-1) - temp_ave;
        //}
        //cout << "iter" << iter << endl;
        //cout << temp_ave << endl;
        //cout << mean(alpha_em) << endl;
        //end test
        
        
        for(int i=1; i<=n_em; i++){
            arma::mat Bi_em     = sub_m(B_em, ni_em, i);
            arma::vec zi_em     = sub_v(z_em, ni_em, i);

            arma::vec alphai_em = alpha_em.col(i-1);
            arma::mat Sigmai_em = Sigma_em.rows(J_em*(i-1),J_em*i-1);
            
            Mi_em = alphai_em * alphai_em.t() + Sigmai_em;
            M_em.rows(J_em*(i-1), J_em*i - 1) = Mi_em;
            temp2 = temp2 + Mi_em;
            arma::mat s2_em = s2_update(Theta_em, thetamu_em, zi_em, Bi_em, alphai_em, Sigmai_em);
            temp = temp + s2_em(0);   //change 1*1 matrix into double:(0)
            temp3 = temp3 + Bi_em.t() * Bi_em;
            temp4 = temp4 + Bi_em.t() * (zi_em - Bi_em * Theta_em * alphai_em);
        }  //update alpha,Sigma, M and thetamu
        
        sigma2hat_em  = temp/sum(ni_em);
        Dhat_em       = temp2/n_em;
        thetamuhat_em = inv(temp3 + sigma2_em * lambmu_em * P_em) * temp4;
        //test Aug 12, 2021
        
        for(int j=1; j<=J_em; j++){
            arma::mat temp5 = zeros(nbasis_em,nbasis_em);
            arma::vec temp7 = zeros(nbasis_em);
            for(int i=1; i<=n_em; i++){
                arma::mat Bi_em     = sub_m(B_em, ni_em, i);
                arma::vec zi_em     = sub_v(z_em, ni_em, i);
                arma::mat Sigmai_em = Sigma_em.rows(J_em*(i-1), J_em*i - 1);
                
                temp5 = temp5 + (alpha_em(j-1,i-1)*alpha_em(j-1,i-1) + Sigmai_em(j-1,j-1)) * Bi_em.t() * Bi_em;
                
                arma::vec temp6 = zeros(ni_em(i-1));
                for(int k=1; (k!=j) && k <= J_em; k++){
                    temp6 = temp6 + (alpha_em(k-1,i-1)*alpha_em(j-1,i-1) + Sigmai_em(k-1,j-1)) * Bi_em * Theta_em.col(k-1);
                }
                
                temp7 = temp7 + Bi_em.t() * ((zi_em - Bi_em*thetamu_em)*alpha_em(j-1,i-1) - temp6);
            }
            Thetahat_em.col(j-1) = inv(temp5 + sigma2_em * lambpc_em * P_em) * temp7;
        }  //update thetaj, or say Theta
        
        
        //testing
        
        arma::mat Thetatilde_em = Theta_em;
            
        eig_sym(eigvalue, eigvec, Thetahat_em * Dhat_em * Thetahat_em.t());
        //since the eigenvalue decomposition result arrange from small to big,
        //need to chance the order, using flip_v and flip_m, functions defined.
        
        Dhat_em = diagmat(flip_v(eigvalue.subvec(nbasis_em - J_em, nbasis_em - 1)));
        Thetahat_em = flip_m(eigvec.cols(nbasis_em - J_em, nbasis_em - 1));
        
 
        arma::mat trans_alpha_em = Thetahat_em.t() * Thetatilde_em;
        
        for(int t=1; t<= n_em; t++){
            alpha_em.col(t-1) = trans_alpha_em * alpha_em.col(t-1);
        }
        
        
        //for(int j = 1; j <= J_em; j++){
        //    if(Thetahat_em(0,j-1)<0) Thetahat_em.col(j-1) = -1 * Thetahat_em.col(j-1);
        //}
        
        sigma2_em  = sigma2hat_em;  //update sigma2
        D_em       = Dhat_em;       //update D
        thetamu_em = thetamuhat_em; //update thetamu
        Theta_em   = Thetahat_em;   //update Theta


    } //iter loop
    
//    Theta_em = - Theta_em;
//    alpha_em = - alpha_em;  //July 21, 2020, delete this line.
    //int iiii = 1;
    return List::create(Named("thetamu")  = thetamu_em,
                        Named("BigTheta") = Theta_em,
                        Named("alpha")    = alpha_em,
                        Named("sigma2")   = sigma2_em,
                        Named("iter")     = iter,
                        Named("D")        = D_em,
                        Named("Sigma")    = Sigma_em);
    
}

//select submatrix of B, etc.
arma::mat sub_m(const arma::mat& A_m, const arma::vec& n_m, const int& i){
    
    arma::mat B_m;
    if(i == 1) B_m = A_m.rows(0, n_m(0)-1);
    else B_m = A_m.rows(sum(n_m.subvec(0,i-2)), sum(n_m.subvec(0,i-1))-1);
    
    return B_m;
}

//select subvector of z, etc.
arma::vec sub_v(const arma::vec& A_v, const arma::vec& n_v, const int& i){
    
    arma::vec B_v;
    if(i == 1) B_v = A_v.rows(0, n_v(0)-1);
    else B_v = A_v.rows(sum(n_v.subvec(0,i-2)), sum(n_v.subvec(0,i-1))-1);
    
    return B_v;
}
                                            

//updating alphai
arma::vec alpha_update(const arma::mat& D_al, const arma::mat& Theta_al, const arma::vec thetamu_al, const double sigma2_al, const arma::mat Bi_al, const arma::vec zi_al){
    
    int ni = zi_al.n_elem;
    arma::vec alphai_al = zeros(D_al.n_rows);
    arma::mat temp = zeros(ni,ni);
    temp = Bi_al * Theta_al * D_al * Theta_al.t() * Bi_al.t() + sigma2_al * eye(ni,ni);
    temp = inv(temp);
    alphai_al = D_al * Theta_al.t() * Bi_al.t() * temp * (zi_al - Bi_al * thetamu_al);
    
    return alphai_al;
}

//updating Sigmai
arma::mat Sigma_update(const arma::mat& D_si, const arma::mat& Theta_si, const double sigma2_si, const arma::mat& Bi_si){
    
    int nr = D_si.n_rows;
    int ni = Bi_si.n_rows;
    arma::mat Sigma_si = zeros(nr,nr);
    arma::mat temp = zeros(ni,ni);
    temp = Bi_si * Theta_si * D_si * Theta_si.t() * Bi_si.t() + sigma2_si * eye(ni,ni);
    temp = inv(temp);
    Sigma_si = D_si - D_si * Theta_si.t() * Bi_si.t() * temp * Bi_si * Theta_si * D_si;
    
    return Sigma_si;
}

//for updating part of sigma2
arma::mat s2_update(const arma::mat& Theta_s, const arma::vec& thetamu_s, const arma::vec& zi_s, const arma::mat& Bi_s, const arma::vec& alphai_s, const arma::mat& Sigmai_s){
    
    arma::mat sig;
    arma::vec temp = zeros(zi_s.n_elem);
    temp = zi_s - Bi_s * thetamu_s - Bi_s * Theta_s * alphai_s;
    sig = temp.t() * temp + trace(Bi_s * Theta_s * Sigmai_s * Theta_s.t() * Bi_s.t());
    
    return sig;
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

