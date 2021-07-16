//[[Rcpp::depends( RcppArmadillo )]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// #include <cmath>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

double T_fun(const double& x, const int& type){
  double act;
  if(type == 1){
    if(x > 0.0){
      act = 1.0;
    }else{
      act = 0.0;
    }
  }
  if(type == 2){
    if(x > 0.0){
      act = x;
    }else{
      act = 0.0;
    }
  }
  return (act);
}


// type 1 : Step function 0 or 1
// type 2 : ReLU; max(t,alpha0)

// [[Rcpp::export]]
Rcpp::List BIVS(const arma::vec & y, const arma::mat & X, int  N, double N1,
                           int BURN, arma::colvec alpha, arma::colvec zeta,arma::colvec w,
                           arma::colvec beta, arma::colvec z, arma::colvec sig,
                           double n1, double p1, double  tau_w, double  tau_z,
                           double  alpha0, double  beta0, double size_a,
                           double size_b, const int & sig_update,const int & zeta_update,
                           const int & K, const int & L, const int& type){
  int n =  X.n_rows, p = X.n_cols;
  arma::mat THETA(n,p), GAM(n,p);
  arma::colvec theta_sj(n);
  arma::colvec u(n);
  arma::colvec theta0(p);
  arma::colvec y_hat(n), y_hat_save(n);
  arma::colvec res(n), rx(1), r(1);
  arma::colvec cand(1), curr(1);
  arma::colvec theta_sj_cand(n), theta_sj_curr(n);
  arma::colvec rx_cand(1), rx_curr(1);
  double a_sq_cand, a_sq_curr;
  arma::colvec beta_cand(p);
  arma::colvec alpha_cand1(p);
  arma::colvec z_cand(p);
  //arma::colvec zeta_cand(p);
  arma::colvec u_cand(u);
  arma::colvec u_curr(u);
  double T_a,  T_b;
  arma::colvec res_cand(n), res_curr(n);
  arma::colvec gam_alpha(p), gam_beta(p);
  double alpha_cand, w_cand, zeta_cand;
  double a_sq;
  double aa, bb;
  double mu1, S;
  double ratio;
  double pq = 0.0;
  double pq_bz = 0.0;
  double obj;
  double beta00;
  IntegerVector sam(p), sam1(p);
  arma::colvec X_norm(p);
  arma::mat SAVE_zeta(p,N); SAVE_zeta.zeros();
  arma::mat SAVE_alpha(p,N); SAVE_alpha.zeros();
  arma::mat SAVE_w(p,N); SAVE_w.zeros();
  arma::mat SAVE_beta(p,N); SAVE_beta.zeros();
  arma::mat SAVE_z(p,N); SAVE_z.zeros();
  arma::colvec OBJ(N), SIG(N);
  arma::mat THETA_save(n,p);
  arma::mat GAM_save(n,p);
  arma::colvec theta2_save(p), theta2(p);
  arma::colvec gam2_save(p), gam2(p);
  //Timer timer;
  int i;
  int i0;
  int ii;
  int j,jj, jj0;
  int k;
  int i1,i2;
  //int start;
  //double pq = 0.0;
  
  for(j=0; j<p; j++){
    sam1(j) = j;
    zeta(j) = 1.0;
    X_norm(j) = sum(X.col(j) % X.col(j));
    if(beta(j)>beta0){
      theta0(j) = (beta(j)-beta0)*z(j);
    }else{
      theta0(j) = 0.0;
    }
    gam_alpha(j) = 0.0;
    gam_beta(j) = 0.0;
  }
  u = X*theta0;
  arma::uvec ind_beta = arma::find(gam_beta == 1.0);
  arma::uvec ind_alpha = arma::find(gam_alpha == 1.0);
  
  for(ii=0; ii<n; ii++){
    for(j=0; j<p; j++){
      aa = alpha(j) + u(ii)*zeta(j)  - alpha0;
      if( aa > 0.0 ){
        THETA(ii,j) =  w(j);
      }else{
        THETA(ii,j) = 0.0;
      }
      THETA_save(ii,j) = 0.0;
      GAM_save(ii,j) = 0.0;
      theta2_save(j) = 0.0;
      gam2_save(j) = 0.0;
    }
  }
  
  
  for(ii=0; ii<n; ii++){
    aa = 0;
    for(j=0; j<p; j++){
      aa = aa + X(ii,j)*THETA(ii,j);
    }
    y_hat(ii) = aa;
  }
  Timer timer;
  //  alpha00 = alpha0;
  beta00 = beta0;
  /////////////
  // start MCMC
  for(i=0; i<(N+BURN); i++){
    beta0 = beta00;
    res = y - y_hat;
      for(ii=0; ii<p; ii++){
        j = sam(ii);
        res = res +  (X.col(j) % THETA.col(j));
        rx(0) = 0.0;
        a_sq = 0.0;
        for(k=0; k<n; k++){
          aa = alpha(j) + u(k)*zeta(j) - alpha0;
          if( aa > 0 ){
            rx(0) = rx(0) + res(k)*X(k,j);
            //theta_sj(k) = 1.0;
            a_sq = a_sq + X(k,j)*X(k,j);
          }else{
            //theta_sj(k) = 0.0;
          }
          theta_sj(k) = T_fun(aa,type);
        }
        S = 1.0 / (a_sq + 1.0/tau_w);
        mu1 = rx(0)*S;
        r(0) = R::rnorm( 0.0 , 1.0 );
        w(j) = mu1 + sqrt(sig(0)*S)*r(0);
        THETA.col(j) = theta_sj*w(j);
        res = res -  (X.col(j) % THETA.col(j));
      }
      y_hat = y - res;
    
    res = y - y_hat;
    // sample alpha and w
    sam = RcppArmadillo::sample(sam1, p, FALSE, NumericVector::create());
    for(ii=0; ii<p; ii++){
      j = sam(ii);
      for(i0=0; i0<K; i0++){
        res = res +  (X.col(j) % THETA.col(j));
        
        r(0) = R::rnorm( 0.0 , 1.0 );
        alpha_cand = alpha(j) + size_a*r(0);
        r(0) = R::rnorm( 0.0 , 1.0 );
        w_cand = w(j);// + size_a*r(0)/10.0;
        
        // evaluate the current status
        rx_curr(0) = 0.0;
        a_sq_curr = 0.0;
        for(k=0; k<n; k++){
          aa = alpha(j) + u(k)*zeta(j) - alpha0;
          if( aa > 0.0 ){
            rx_curr(0) = rx_curr(0) + res(k)*X(k,j)*w(j);
            theta_sj_curr(k) = T_fun(aa,type)*w(j);
            a_sq_curr = a_sq_curr + X(k,j)*X(k,j)*w(j)*w(j);
          }else{
            //theta_sj_curr(k) = 0.0;
            theta_sj_curr(k) = T_fun(aa,type)*w(j);
          }
          
        }
        //curr = -0.5*log(a_sq_curr + 1.0/tau_w) - 0.5*alpha(j)*alpha(j) + 0.5*rx_curr*rx_curr/(a_sq_curr + 1.0/tau_w)/sig(0);
        curr(0) = -0.5*(a_sq_curr - 2*rx_curr(0))/sig(0) - 0.5*alpha(j)*alpha(j) - 0.5*w(j)*w(j)/(sig(0)*tau_w);
        
        rx_cand(0) = 0.0;
        a_sq_cand = 0.0;
        for(k=0; k<n; k++){
          aa = alpha_cand + u(k)*zeta(j) - alpha0;
          if( aa > 0.0 ){
            rx_cand(0) = rx_cand(0) + res(k)*X(k,j)*w_cand;
            //theta_sj_cand(k) = w_cand;
            theta_sj_cand(k) = T_fun(aa,type)*w_cand;
            a_sq_cand = a_sq_cand + X(k,j)*X(k,j)*w_cand*w_cand;
          }else{
            theta_sj_cand(k) = T_fun(aa,type)*w_cand;
          }
        }
        //curr = -0.5*log(a_sq_curr + 1.0/tau_w) - 0.5*alpha(j)*alpha(j) + 0.5*rx_curr*rx_curr/(a_sq_curr + 1.0/tau_w)/sig(0);
        cand(0) = -0.5*(a_sq_cand - 2*rx_cand(0))/sig(0) - 0.5*alpha_cand*alpha_cand - 0.5*w_cand*w_cand/(sig(0)*tau_w);
        
        aa = R::runif( 0.0, 1.0);
        ratio = cand(0) - curr(0);
        if(ratio > log(aa) ){
          alpha(j) = alpha_cand;
          w(j) = w_cand;
          THETA.col(j) = theta_sj_cand;
          pq = pq + 1.0;
        }
        res = res -   (X.col(j) % THETA.col(j));
      }
      
      y_hat = y - res;
    }
    for(ii=0; ii<n; ii++){
      aa = 0.0;
      for(j=0; j<p; j++){
        aa = aa + X(ii,j)*THETA(ii,j);
      }
      y_hat(ii) = aa;
    }
    res = y - y_hat;
    
    // sample beta and z
    u_cand = u;
    sam = RcppArmadillo::sample(sam1, p, FALSE, NumericVector::create());
    for(i0=0; i0<L; i0++){
      beta_cand = beta;
      z_cand = z;
      for(jj0=0; jj0<p; jj0++){
        jj = sam(jj0);
        
        beta_cand(jj) = beta(jj) + R::rnorm( 0.0 , 1.0 )*size_b;
        z_cand(jj) = z(jj) + R::rnorm( 0.0 , 1.0 )*size_b/10.0;
        // update u
        for(ii=0; ii<n; ii++){
          u_cand(ii) = 0;
          for(j=0; j<p; j++){
            if(beta_cand(j) > beta0){
              T_b = beta_cand(j) - beta0;
              u_cand(ii) = u_cand(ii) + X(ii,j)*T_b*z_cand(j);
            }else{
              T_b = 0.0;
            }
          }
        }
        for(ii=0; ii<n; ii++){
          aa = 0;
          for(j=0; j<p; j++){
            aa = alpha(j) + u_cand(ii)*zeta(j) - alpha0;
            if(aa > 0.0){
              T_a =  T_fun(aa,type)*w(j);
              aa = aa + X(ii,j)*T_a*w(j);
            }else{
              T_a =  T_fun(aa,type)*w(j);
            }
          }
          res_cand(ii) = y(ii) - aa;
        }
        cand(0) = -0.5*sum(res_cand % res_cand)/sig(0) - 0.5*beta_cand(jj)*beta_cand(jj) - 0.5*z_cand(jj)*z_cand(jj)/(tau_z*sig(0));
        res_cand = y - y_hat;
        curr(0) = -0.5*sum(res_cand % res_cand)/sig(0) - 0.5*beta(jj)*beta(jj) - 0.5*z(jj)*z(jj)/(tau_z*sig(0));
        
        // MH step
        aa = R::runif( 0.0, 1.0);
        ratio = cand(0) - curr(0);
        if(ratio > log(aa) ){
          pq_bz = pq_bz + 1.0;
          beta(jj) = beta_cand(jj);
          z(jj) = z_cand(jj);
          u = u_cand;
          // update THETA
          for(ii=0; ii<n; ii++){
            for(j=0; j<p; j++){
              aa = alpha(j) + u(ii)*zeta(j) - alpha0;
              if( aa > 0.0 ){
                THETA(ii,j) =  w(j);
              }else{
                THETA(ii,j) = 0.0;
              }
            }
          }
          // update y_hat and residuals
          for(ii=0; ii<n; ii++){
            aa = 0.0;
            for(j=0; j<p; j++){
              aa = aa + X(ii,j)*THETA(ii,j);
            }
            y_hat(ii) = aa;
          }
          res = y - y_hat;
        }else{
          beta_cand(jj) = beta(jj);
          z_cand(jj) = z(jj);
        }
      }
    }
    res = y - y_hat;
    // update zeta

    for(j=0; j<p; j++){
      if(beta(j) > beta0){
        T_b = beta(j) - beta0;
        gam2(j) = 1.0;
      }else{
        T_b = 0.0;
        gam2(j) = 0.0;
      }
      theta2(j) =  T_b*z(j);
    }
    for(ii=0; ii<n; ii++){
      aa = 0.0;
      for(j=0; j<p; j++){
        aa = aa + X(ii,j)*THETA(ii,j);
      }
      y_hat(ii) = aa;
    }
    res = y - y_hat;
    
    //update sig(0)
    if( sig_update == 1 ){
      if(i % 10 == 9){
        aa = (n1+p1+p1)*0.5 + 1.0;
        bb = sum(res % res)*0.5 +  sum(w % w)*0.5/tau_w + sum(z % z)*0.5/tau_z;
        sig(0) = R::rgamma(aa,1.0/bb);
        sig(0) = 1.0/sig(0);
      }
    }
    for(ii=0; ii<n; ii++){
      for(j=0; j<p; j++){
        aa = alpha(j) + u(ii) - alpha0;
        if( aa > 0 ){
          GAM(ii,j) =  1.0;
        }else{
          GAM(ii,j) = 0.0;
        }
      }
    }
    
    // update GAM
    obj = -0.5*sum(res % res)/sig(0) -  sum(alpha % alpha)*0.5 -  sum(w % w)*0.5/(tau_w*sig(0)) - 0.5*(n1+p1+2+p1)*log(sig(0));
    obj = obj - sum(beta % beta)*0.5 -  sum(z % z)*0.5/(tau_z*sig(0));
    if((i+1) > BURN){
      OBJ(i-BURN) = obj;
      timer.step("");
      y_hat_save = y_hat_save + y_hat/N1;
      THETA_save = THETA_save + THETA/N1;
      GAM_save = GAM_save + GAM/N1;
      theta2_save = theta2_save + theta2/N1;
      gam2_save = gam2_save + gam2/N1;
      SAVE_alpha.col(i-BURN) = alpha;
      SAVE_w.col(i-BURN) = w;
      SAVE_beta.col(i-BURN) = beta;
      SAVE_z.col(i-BURN) = z;
      SAVE_zeta.col(i-BURN) = zeta;
      SIG(i-BURN) = sig(0);
    }
    if((i+1) % 10 == 0){
      for(j=0; j<p; j++){
        aa = 0;
        for(ii=0; ii<n; ii++){
          if(GAM(ii,j) == 1.0){
            aa = aa + 1.0/n1;
          }
        }
        if(aa > 0.05){
          gam_alpha(j) = 1.0;
        }else{
          gam_alpha(j) = 0.0;
        }
        if(beta(j) > beta0){
          gam_beta(j) = 1.0;
        }else{
          gam_beta(j) = 0.0;
        }
      }
      ind_beta = arma::find(gam_beta == 1.0);
      ind_alpha = arma::find(gam_alpha == 1.0);
      Rcpp::Rcout << "########### RWMH ###########" << std::endl;
      Rcpp::Rcout << "#iterations = " << i+1 << std::endl;
      Rcpp::Rcout << "log-posterior = " << obj << std::endl;
      Rcpp::Rcout << "Acceptance rate for alpha = " << pq/(i*p) << std::endl;
      Rcpp::Rcout << "Acceptance rate for beta and z = " << pq_bz/(i*p) << std::endl;
      Rcpp::Rcout << "Selected variables for global level = " << (ind_alpha+1).t() << std::endl;
      Rcpp::Rcout << "Selected variables for local level = " << (ind_beta+1).t() << std::endl;
      Rcpp::Rcout << "alpha0 = " << alpha0 << std::endl;
      Rcpp::Rcout << "beta0 = " << beta0 << std::endl;
      Rcpp::Rcout << "sigma^2 = " << sig(0) << std::endl;
    }
    
    
  }
  return Rcpp::List::create(Rcpp::Named("y_hat") = y_hat_save,
                            Rcpp::Named("POST") = OBJ,
                            Rcpp::Named("THETA") = THETA_save,
                            Rcpp::Named("theta2") = theta2_save,
                            Rcpp::Named("gam2") = gam2_save,
                            Rcpp::Named("GAM") = GAM_save,
                            Rcpp::Named("alpha") = SAVE_alpha,
                            Rcpp::Named("zeta") = SAVE_zeta,
                            Rcpp::Named("w") = SAVE_w,
                            Rcpp::Named("beta") = SAVE_beta,
                            Rcpp::Named("z") = SAVE_z,
                            Rcpp::Named("sig") = SIG,
                            Rcpp::Named("u") = u,
                            Rcpp::Named("CompTime") = timer
                              //Rcpp::Named("accpt") = pq
  );
  
}
