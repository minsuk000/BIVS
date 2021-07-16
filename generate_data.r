# this is the file for the code that generates data from well specified BIVS

T_step = function(x){
  res = c()
  for(i in c(x)){
    if(i>0){
      res = c(res, 1)
    } else{
      res = c(res, 0)
    }
  }
  return(res)
}

T_ReLU = function(x){
  res = c()
  for(i in c(x)){
    if(i>0){
      res = c(res, i)
    } else{
      res = c(res, 0)
    }
  }
  return(res)
}


generate_THETA_star = function(bet, bet0, z, T_function){
  if(T_function == 'step')
    theta_star = T_step(bet-bet0)*z
  if(T_function == 'relu')
    theta_star = T_ReLU(bet-bet0)*z
  return(theta_star)
}

generate_THETA = function(alp, alp0, zeta, w, theta_star, covariates, T_function){
    uu = covariates%*%theta_star
    # KK = T0(tcrossprod(rep(1,dim(covariates)[1]),alp) + tcrossprod(uu,zeta) - alp0)
    THETAA = covariates
    for(ii in 1:(dim(covariates)[1])){
      for(jj in 1:(dim(covariates)[2])){
        if(T_function == 'step'){
          THETAA[ii,jj] = T_step(alp[jj]+uu[ii]*zeta[jj]-alp0)*w[jj]
        }
        if(T_function == 'relu'){
          THETAA[ii,jj] = T_ReLU(alp[jj]+uu[ii]*zeta[jj]-alp0)*w[jj]
        }
      }
    }
    return(THETAA)
}
  
generate_mean = function(alp, alp0, zeta, w, theta_star, covariates, T_function){
    THETAA = generate_THETA(alp, alp0, zeta, w, theta_star, covariates, T_function)
    return(apply(covariates*THETAA,1,sum))
}
