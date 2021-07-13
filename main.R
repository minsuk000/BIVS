# main file

library(NPrior)
library(hierNet)
library(glmnet)
library(MASS)
source(generate_data.r)

BURN = 2000
N = 5000
n = 500
p = 10
signal_level = 3
T_function = "step" # choices are "step" and "relu"
data_pre_process = FALSE
set.seed(91817)
nsplit <- 16 # number of splits to test for MSE


BIVS_boston = function(i0, tr, xall, yall, X1_intr, N, BURN, i){
  X0 = xall
  y0 = yall
  n = nrow(xall)
  # print(n)
  otr <- tr[[i0]] # training set
  ind_test = (1:n)[-otr] # test set
  n_test = length(ind_test)
  # print(n_test)
  n = n - n_test # n is not |training set|
  # print(n)
  X = X0[-ind_test,]
  y = matrix(y0[-ind_test],n,1)
  X_test = X0[ind_test,]
  y_test = matrix(y0[ind_test],n_test,1)
  X_intr = X1_intr[-ind_test,]
  X_intr_test = X1_intr[ind_test,]
  
  
  xtr = X
  ytr = y
  xte = X_test
  yte = y_test
  
  n = nrow(X)
  p = ncol(X)
  
  # try fully specified first
  # a0_cand = c(0.5,1,1.3,1.5,1.7,2,2.2,2.5,3,3.4,3.7,4)
  # a0_cand = c(-a0_cand, a0_cand) # values in CV
  # try fully specified first
  
  a0_cand = 0
  sig = matrix(1,1,1)
  
  K = 5 # CV fold to pick parameters
  rspe_bivs0 = post = rep(0,5)
  size_bivs0 = rep(0,5)
  # print(n)
  size = floor(n/K)
  # print(size)
  cv.err = matrix(0,K,length(a0_cand))
  for(v in 1:length(a0_cand)){
    a0 = b0 = a0_cand[v]
    for(h in 1:K){
      ind_val = ((h-1)*size+1):(h*size) # validation index
      # print(ind_val)
      X_tr = X[-ind_val,]
      X_val = X[ind_val,]
      y_tr = y[-ind_val,]
      y_val = y[ind_val,]
      beta = rnorm(p)
      alpha = rnorm(p)
      w = rnorm(p)
      z = rnorm(p)
      zeta= rnorm(p)
      beta = rep(0,p);#rnorm(p)
      alpha = rep(0,p);#rnorm(p)
      w = rep(0,p);#rnorm(p)
      z = rep(0,p);#rnorm(p)
      # zeta = rep(1,p);#rep(1,p)          ##??
      zeta = rnorm(p)
      sig = matrix(1,1,1)
      
      pmt = proc.time()[3]
      
      # fit.RW = IVS_RW_naive000(y_tr, X_tr,  N, N1=N,
      #                          BURN, alpha_true, zeta_true,
      #                          w_true, beta_true, z_true,
      #                          sig, n1=n, p1=p, tau_w=1,
      #                          tau_z=1, alpha0=a0, beta0=b0,
      #                          size_a=3, size_b=3, sig_update=1,
      #                          K = 10, L = 10, zeta_update = 1)
      # WE SHOULD USE JUST ONE CPP FUNCTION FOR FITTING, SPECIFYING OTHERS AS ARGUMENTS
      fit.RW = BIVS('RW', 'naive', y_tr, X_tr,  N, N1=N,
                    BURN, alpha_true, zeta_true,
                    w_true, beta_true, z_true,
                    sig, n1=n, p1=p, tau_w=1,
                    tau_z=1, alpha0=a0, beta0=b0,
                    size_a=3, size_b=3, sig_update=1,
                    K = 10, L = 10, zeta_update = 1)
      pmt0 = proc.time()[3] - pmt
      # THETA = fit.RW$THETA
      # print(X_val)
      # n_val = nrow(X_val)
      # print(n_val)
      # print(p)
      # theta_test = matrix(0, n_val,p)
      # alpha = apply(fit.RW$alpha,1,mean)
      # beta = apply(fit.RW$beta,1,mean)
      # w =  apply(fit.RW$w,1,mean)
      # z =  apply(fit.RW$z,1,mean)
      # beta1 = ifelse(beta-b0>0,beta-b0,0)*z
      # u_test = X_val%*%beta1
      # for(j in 1:p){
      #   thta = ifelse(u_test + alpha[j] > a0,1,0)*w[j]
      #   theta_test[,j] = thta
      # }
      # ind = which(fit.RW$GAM>0.5,arr.ind=T)
      # ind = unique(ind[,2])
      # size_bivs = length(ind)
      # pred_bivs = apply(X_val*theta_test,1,sum) # is this the right way to test?
      
      pred_bivs = BIVS_posterior_mean(fit.RW, X_val) # THIS FUNCTION HAS NOT BEEN WRITTEN YET
      cv.err[h,v] = sqrt(mean((pred_bivs-y_val)^2))
    }
    # print(v)
  }
  cv.err0 = apply(cv.err,2,mean)
  a0 = b0 = a0_cand[ which.min(cv.err0) ]
  
  
  plot(cv.err0,type="l")
  
  # Once a0, b0 are chosen, use multiple starts to fit
  # WE ARE NOW DOING WELL SPECIFIED, SO THERE IS NO DIFFERENCE, BUT IN THE FUTURE, SHOULD USE MULTIPLE STARTS IN 
  # THE CV STEP AS WELL
  for(h in 1:5){ # doing 5 multiple starts
    # CAN WE INTEGRATE THE MULTIPLE START FUNCTIONALITY IN THE CPP FUNCTION?
    beta = rnorm(p)
    alpha = rnorm(p)
    w = rnorm(p)
    z = rnorm(p)
    zeta= rnorm(p)
    beta = rep(0,p);#rnorm(p)
    alpha = rep(0,p);#rnorm(p)
    w = rep(0,p);#rnorm(p)
    z = rep(0,p);#rnorm(p)
    zeta= rep(0,p);#rep(1,p)
    sig = matrix(1,1,1)
    # fit.RW = IVS_RW_new(y, X,y,X,N,N1=N,BURN,alpha,w,beta, z,gamma,v,sig, n1=n, p1=p, tau_w=1,tau_z=1,tau_v=1,
    #                                        alpha0, beta0, gamma0,size_a=3, size_b=3, sig_update=1, zeta_update=0, K=5, L=5)
    fit.RW = BIVS('RW', 'naive', y, X, N, N1=N,
                             BURN, alpha, zeta,
                             w, beta, z,
                             sig, n1=n, p1=p, tau_w=1,
                             tau_z=1, alpha0=a0, beta0=b0,
                             size_a=3, size_b=3, sig_update=1,
                             K = 10, L = 10, zeta_update = 1)
    THETA = fit.RW$THETA
    theta_test = matrix(0, n_test, p)
    alpha = apply(fit.RW$alpha,1,mean)
    beta = apply(fit.RW$beta,1,mean)
    w =  apply(fit.RW$w,1,mean)
    z =  apply(fit.RW$z,1,mean)
    beta1 = ifelse(beta-b0>0,beta-b0,0)*z
    u_test = X_test%*%beta1
    for(j in 1:p){
      thta = ifelse(u_test + alpha[j] > a0,1,0)*w[j]
      theta_test[,j] = thta
    }
    ind = which(fit.RW$GAM>0.5,arr.ind=T)
    ind = unique(ind[,2])
    size_bivs0[h] = length(ind)
    theta = apply(X_test*theta_test,1,sum)
    rspe_bivs0[h] = sqrt(mean((theta-y_test)^2))
    post[h] = mean(fit.RW$POST) 
  }
  ind = which.max(post)
  size_bivs1 = size_bivs0[ind]
  rspe_bivs1 = rspe_bivs0[ind]
  ##################################################
  
  size = floor(n/K)
  cv.err = matrix(0,K,length(a0_cand))
  for(v in 1:length(a0_cand)){
    a0 = b0 = a0_cand[v]
    for(h in 1:K){
      ind_val = ((h-1)*size+1):(h*size)
      X_tr = X[-ind_val,]
      X_val = X[ind_val,]
      y_tr = y[-ind_val,]
      y_val = y[ind_val,]
      beta = rnorm(p)
      alpha = rnorm(p)
      w = rnorm(p)
      z = rnorm(p)
      zeta= rnorm(p)
      beta = rep(0,p);#rnorm(p)
      alpha = rep(0,p);#rnorm(p)
      w = rep(0,p);#rnorm(p)
      z = rep(0,p);#rnorm(p)
      zeta= rep(0,p);#rep(1,p)
      sig = matrix(1,1,1)
      
      pmt = proc.time()[3]
      fit.RW = BIVS('RW', 'simplest', y_tr, X_tr, N ,N1=N,BURN,alpha,zeta,w,beta, z,sig, n1=n, p1=p, tau_w=1,tau_z=1,
                               a0, b0, size_a=3, size_b=3, sig_update=1, zeta_update=0,
                               K=10, L=10)
      pmt0 = proc.time()[3] - pmt
      THETA = fit.RW$THETA
      n_val = nrow(X_val)
      theta_test = matrix(0, n_val,p)
      alpha = apply(fit.RW$alpha,1,mean)
      beta = apply(fit.RW$beta,1,mean)
      w =  apply(fit.RW$w,1,mean)
      z =  apply(fit.RW$z,1,mean)
      beta1 = ifelse(beta-b0>0,beta-b0,0)*z
      u_test = X_val%*%beta1
      for(j in 1:p){
        thta = ifelse(u_test + alpha[j] > a0,1,0)*w[j]
        theta_test[,j] = thta
      }
      ind = which(fit.RW$GAM>0.5,arr.ind=T)
      ind = unique(ind[,2])
      size_bivs = length(ind)
      pred_bivs = apply(X_val*theta_test,1,sum)
      cv.err[h,v] = sqrt(mean((pred_bivs-y_val)^2))
    }
    print(v)
  }
  cv.err0 = apply(cv.err,2,mean)
  a0 = b0 = a0_cand[ which.min(cv.err0) ]
  plot(cv.err0,type="l")
  
  for(h in 1:5){
    beta = rnorm(p)
    alpha = rnorm(p)
    w = rnorm(p)
    z = rnorm(p)
    zeta= rnorm(p)
    beta = rep(0,p);#rnorm(p)
    alpha = rep(0,p);#rnorm(p)
    w = rep(0,p);#rnorm(p)
    z = rep(0,p);#rnorm(p)
    zeta= rep(0,p);#rep(1,p)
    sig = matrix(1,1,1)
    fit.RW = BIVS('RW', 'simplest', y, X,  N, N1=N, BURN, alpha, zeta, w,
                             beta,  z, sig, n1 = n, p1 = p, tau_w=1, tau_z=1,
                             alpha0=a0, beta0=b0, size_a=3, size_b=3, sig_update=1, zeta_update=0,
                             K=10, L=10)
    THETA = fit.RW$THETA
    theta_test = matrix(0, n_test, p)
    alpha = apply(fit.RW$alpha,1,mean)
    beta = apply(fit.RW$beta,1,mean)
    w =  apply(fit.RW$w,1,mean)
    z =  apply(fit.RW$z,1,mean)
    beta1 = ifelse(beta-b0>0,beta-b0,0)*z
    u_test = X_test%*%beta1
    for(j in 1:p){
      thta = ifelse(u_test + alpha[j] > a0,1,0)*w[j]
      theta_test[,j] = thta
    }
    ind = which(fit.RW$GAM>0.5,arr.ind=T)
    ind = unique(ind[,2])
    size_bivs0[h] = length(ind)
    theta = apply(X_test*theta_test,1,sum)
    rspe_bivs0[h] = sqrt(mean((theta-y_test)^2))
    post[h] = mean(fit.RW$POST) 
  }
  ind = which.max(post)
  size_bivs = size_bivs0[ind]
  rspe_bivs = rspe_bivs0[ind]
  ##################################################
  #  IVS_interact_RW(y, X,  N, N1=N, BURN, alpha, zeta, w,
  #                  beta,  z, sig, n1 = n, p1 = p, tau_w=1, tau_z=1,
  #                  alpha0=a0, beta0=b0, size_a=3, size_b=3, sig_update=1, zeta_update=0,
  #                  K=10, L=10)
  print('a')
  fit = NPrior_run(X, y, alpha0_update = F, alpha0= a0, sig_update=T, sig=1)
  beta_lm = fit$ThetaHat # apply(fit$THETA,1,mean)
  ind = which(apply(fit$AlphaSamples,1,mean) > a0)
  size_lm = length(ind)
  pred_lm = X_test%*%beta_lm
  ##################################################
  fit = NPrior_run(X_intr, y, alpha0_update = F, alpha0= a0, sig_update=T, sig=1)
  beta_lm_intr = apply(fit$ThetaSamples,1,mean)
  ind = which(apply(fit$AlphaSamples,1,mean) > a0)
  size_lm_intr = length(ind)
  pred_lm_intr = X_intr_test%*%beta_lm_intr
  ##################################################
  p = ncol(X_intr)
  beta = rep(0,p);#rnorm(p)
  alpha = rep(0,p);#rnorm(p)
  w = rep(0,p);#rnorm(p)
  z = rep(0,p);#rnorm(p)
  zeta= rep(0,p);#rep(1,p)
  sig = matrix(1,1,1)
  pmt = proc.time()[3]
  #fit.RW = IVS_RW_simplest(y, X_intr,N,N1=N,BURN,alpha,zeta,w,beta, z,sig, n1=n, p1=ncol(X_intr), tau_w=1,tau_z=1,
  #                         a0, b0, size_a=3, size_b=3, sig_update=1, zeta_update=0,
  #                         K=10, L=10)
  #pmt0 = proc.time()[3] - pmt
  #THETA = fit.RW$THETA
  #theta_test = matrix(0, n_test,p)
  #alpha = apply(fit.RW$alpha,1,mean)
  #beta = apply(fit.RW$beta,1,mean)
  #w =  apply(fit.RW$w,1,mean)
  #z =  apply(fit.RW$z,1,mean)
  #beta1 = ifelse(beta-b0>0,beta-b0,0)*z
  #u_test = X_intr_test%*%beta1
  #for(j in 1:p){
  #  thta = ifelse(u_test + alpha[j] > a0,1,0)*w[j]
  #  theta_test[,j] = thta
  #}
  #ind = which(fit.RW$GAM>0.5,arr.ind=T)
  #ind = unique(ind[,2])
  #size_bivs_intr = length(ind)
  #pred_bivs_intr = apply(X_intr_test*theta_test,1,sum)
  size_bivs_intr = 0 
  ##################################################
  
  zztr <- compute.interactions.c(xtr, diagonal=F)
  zzte <- compute.interactions.c(xte, diagonal=F)
  mzz <- colMeans(zztr)
  zztr <- scale(zztr, center=mzz, scale=FALSE)
  zzte <- scale(zzte, center=mzz, scale=FALSE)
  ztr <- cbind(xtr, zztr)
  zte <- cbind(xte, zzte)
  
  #lam.mel <- exp(seq(log(0.98), log(0.1e-7), length=40))
  #lam.apl <- exp(seq(log(1.36), log(0.05), length=40))
  #mel <- glmnet(xtr, ytr, standardize=FALSE, lambda=lam.mel)
  #phat.mel <- predict(mel, type="response", newx=xte)
  #apl <- glmnet(ztr, ytr, standardize=FALSE, lambda=lam.apl)
  #phat.apl <- predict(apl, type="response", newx=zte)
  #hl <- hierNet.path(xtr, ytr, zz=zztr, minlam=1, maxlam=200, nlam=40, diagonal=F,
  #                   strong=FALSE, maxiter=1000, step=10, trace=2)
  #phat.hl <- predict(hl, newx=xte, newzz=zzte)
  zztr <- compute.interactions.c(X, diagonal=F)
  zzte <- compute.interactions.c(X_test, diagonal=F)
  hl <- hierNet.path(X, as.vector(y), zz=zztr, minlam=1, maxlam=200, nlam=10, diagonal=F,
                     strong=F, maxiter=1000, step=10, trace=2, stand.main=F)
  fitcv=hierNet.cv(hl,X,as.vector(y))
  lamhat = fitcv$lamhat
  fit2=hierNet(X,as.vector(y),lam=lamhat, stand.main=F)
  yhat_hLASSO=predict(fit2,X_test)
  
  ###################################################
  fit_cv = cv.glmnet(X_intr,y)
  lam_hat = fit_cv$lambda.min
  fit_LASSO = glmnet(X_intr,y,lambda = lam_hat)
  beta_LASSO_intr = fit_LASSO$beta
  yhat_LASSO_intr = X_intr_test%*%beta_LASSO_intr
  print(sum((yhat_LASSO_intr - y_test)^2)/n_test)
  ##################################################
  library(glmnet)
  fit_cv = cv.glmnet(X,y)
  lam_hat = fit_cv$lambda.min
  fit_LASSO = glmnet(X,y,lambda = lam_hat)
  beta_LASSO = fit_LASSO$beta
  yhat_LASSO = X_test%*%beta_LASSO
  print(sum((yhat_LASSO - y_test)^2)/n_test)
  
  ##################################################
  #rspe_bivs = sqrt(mean((pred_bivs-y_test)^2))
  #rspe_bivs1 = sqrt(mean((pred_bivs1-y_test)^2))
  
  size_bivs_intr = 0
  rspe_bivs_intr = 0#sqrt(mean((pred_bivs_intr-y_test)^2))
  
  rspe_lm = sqrt(mean((pred_lm-y_test)^2))
  rspe_lm_intr = sqrt(mean((pred_lm_intr-y_test)^2))
  
  return( list(bivs = c(rspe_bivs, size_bivs), bivs1 = c(rspe_bivs1, size_bivs1), bivs_intr = c(rspe_bivs_intr, size_bivs_intr), lm = c(rspe_lm, size_lm) ,lm_intr = c(rspe_lm_intr, size_lm_intr),
               phat.mel=yhat_LASSO, phat.apl=yhat_LASSO_intr,
               phat.hl=yhat_hLASSO, yte= yte) )
}



for(iii in 1:1){
  beta_true = rnorm(p)
  alpha_true = rnorm(p)
  gamma_true = rnorm(p)
  alpha0 = beta0 = 0
  gamma0 = 0
  sig.w = 1
  sig0 = 1
  ones.n = rep(1,n)
  ones.p = rep(1,p)
  w_true = rnorm(p)*sqrt(sig.w*sig0)
  z_true = rnorm(p)*sqrt(sig.w*sig0)
  theta_star_true = generate_THETA_star(beta_true, beta0, z_true, T_function)
  zeta_true = rnorm(p)
  
  X0 = matrix(rnorm(n*p), ncol = p)
  y0 = generate_mean(alpha_true, alpha0, zeta_true, w_true, theta_star_true, X0) + rnorm(n)*sqrt(sig0)
  
  if(data_pre_process){
    X0 = scale(X0)
    y0 = y0-mean(y0)
  }


  # getting all the interactions
  X1_intr = NULL
  for(j in 1:(length(1:p)-1)){
    for(k in (j+1):(length(1:p))){
      X1_intr = cbind(X1_intr, X0[,j]*X0[,k])
      s = ncol(X1_intr)
      colnames(X1_intr)[s] = paste(j,"+",k,sep="")
    }
  }
  
  X1_intr = cbind(X0,X1_intr)
  
  tr <- list()
  for (i in seq(nsplit))
    tr[[i]] <- sample(1:n, size=trunc(n*0.8))
  
  
  # setwd("~/Dropbox/Paper work/IndVarSel/2020/real/hierLASSO/boston")
  # Rcpp::sourceCpp('IVS_RW_simplest.cpp')
  # Rcpp::sourceCpp('IVS_RW_naive000.cpp')
  # BIVS_boston(1, tr, X0, y0, X1_intr,N, BURN, 1)
  library(snowfall)
  sfInit(cpus=16, parallel=T)
  sfLibrary(NPrior)
  sfLibrary(hierNet)
  sfLibrary(glmnet)
  sfLibrary(Rcpp)
  sfExportAll()
  sfClusterEval(Rcpp::sourceCpp('IVS_RW_simplest.cpp'))
  sfClusterEval(Rcpp::sourceCpp('IVS_RW_naive000.cpp'))
  # sfClusterEval(Rcpp::sourceCpp("IVS_RW_new.cpp"))
  # sfClusterEval(Rcpp::sourceCpp("IVS_RW_multiple_N.cpp"))
  
  #sfExportAll()
  #res_ind = 13:18
  #sss = BIVS_boston(1, tr, X0, y0, X1_intr,N, BURN, 1)
  result_boston = sfLapply(seq(nsplit), BIVS_boston, tr, X0, y0, X1_intr,N, BURN, 1)
  save(result_boston, file=sprintf("X_%d_%d_%d_%d.RData", n, p, signal_level, iii))
  sfStop()


  err.hl = err.apl = err.mel  = NULL
  err.bivs = err.bivs1 = err.bivs_intr = err.lm = err.lm_intr = NULL
  rmse <- function(yhat, yte) {
    sqrt(mean((yhat - yte)^2))
    #sqrt(((yhat - yte)^2))
  }
  for (i in seq(nsplit)){
    yte = result_boston[[i]]$yte
    err.hl[i] <- rmse(result_boston[[i]]$phat.hl, yte)
    err.mel[i] <- rmse(result_boston[[i]]$phat.mel, yte)
    err.apl[i] <- rmse(result_boston[[i]]$phat.apl, yte)
    err.bivs[i] <- result_boston[[i]]$bivs[1]
    err.bivs1[i] <- result_boston[[i]][[2]][1]

    err.bivs_intr[i] <- result_boston[[i]]$bivs_intr[1]
    err.lm[i] <- result_boston[[i]]$lm[1]
    err.lm_intr[i] <- result_boston[[i]]$lm_intr[1]
  }

  names = c("LASSO-MEL","LASSO-APL","HLASSO","BIVS","BIVS-ReLU","BVS","BVS+APL")
  png(filename=sprintf("results/X_%d_%d_%d_%d.png", n, p, signal_level, iii),width=1500,height=1000)
  par(mai=c(0.8,0.8,0.1,0.1))
  boxplot(err.mel,err.apl,err.hl, err.bivs, err.bivs1, err.lm, err.lm_intr,boxwex=0.4, names=names, ylab="RMSPE" )
  dev.off()
}

