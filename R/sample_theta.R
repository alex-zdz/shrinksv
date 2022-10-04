sample_theta <- function(mu, phi, sigma, h0, h,  prior_spec ,expert){

  if(expert$block=="1block"){
  
  shrinksv <- prior_spec$shrinksv
  # calculate bT and BT:
  X <- matrix(c(rep(1,Ttot),c(h0,head(h,-1) )), Ttot , 2)
  gamma_varinv <- 1/prior_spec$gamma_var
  phi_varinv <- 1/prior_spec$phi_var
  B0_inv <- diag(c(gamma_varinv,phi_varinv))
  # BT_inv_11 <- Ttot-1 + B0_inv[1,1]
  # BT_inv_12 <- sum(X[,2])
  # BT_inv_22 <- sum(X[,2]*X[,2])
  # BT <- 1/(BT_inv_11*BT_inv_22 - BT_inv_12^2) * matrix(c(BT_inv_22,-BT_inv_12,-BT_inv_12,BT_inv_11), 2,2, byrow = TRUE)

  # old:
  BT <- solve(crossprod(X)+B0_inv)
  # sum2 <- sum(h)
  # sum3 <- c(h0,head(h,-1) )%*%h
  # bT1 <- BT[1,1]*sum2+BT[1,2]*sum3
  # bT2 <- BT[1,2]*sum2+BT[2,2]*sum3
  # bT<-c(bT1,bT2)
  
  bT <- BT%*%t(X)%*%h
  # propose sigma
  cT <- (Ttot-1)/2
  CT <- .5*(sum(h^2) - t(bT)%*%t(X)%*%h)
  #CT <- .5*(sum(c(h0,h)^2) - t(bT)%*%t(X)%*%h)    #old and wrong
  sigma_prop <- sqrt(1/rgamma(1,cT, CT)) # shape = cT , rate = CT parametrization contrary to calling rgamma from c, then it is shape, scale paramet.

  # Calculate Cholesky decomposition of sig2BT- covariance matrix of gamma and phi:
  chol11 <- sqrt(BT[1,1])*sigma
  chol12 <- (BT[1,2]/chol11)*sigma
  chol22 <- sqrt(BT[2,2]-chol12^2)*sigma

  # propose phi using inverse transform sampling from a truncated normal distribution:
  quant_low <-  pnorm(-1, bT[2], chol22, lower.tail = TRUE ,log.p = FALSE)  # use chol11 = sqrt(BT[1,1]*sigma^2)
  quant_high <- pnorm(1, bT[2], chol22,  lower.tail = TRUE ,log.p = FALSE)
  phi_prop <- qnorm((quant_low + runif(1)*(quant_high-quant_low)), bT[2],sd= chol22,lower.tail = TRUE ,log.p = FALSE)

  # propose gamma
  gamma_prop <- rnorm(1, bT[1] + chol12*((phi_prop-bT[2])/chol22), sd = chol11) # rewrite in terms of BT?

  #03.08.2022 try out different sampling method for phi
  # phi_prop <- bT[2] + chol22*rnorm(1)  #+0.01
  # if (abs(phi_prop)>=1){  #outside the unit sphere
  # phi_prop <- Inf
  # }
  # works but the same problem of slightly underestimating the true value for phi ~ 0.9
  
  
  # Calculate acceptance rates
  # Start with gamma and phi
  mu_prop <- gamma_prop/(1-phi_prop)
  #phi_const <- 1-phi
  proposal_gamma_sd <- sqrt(prior_spec$gamma_var)
  proposal_phi_sd <- sqrt(prior_spec$phi_var)

  # first the likelihood
  ar_prob <- dnorm(h0, mu_prop, sigma_prop/sqrt(1-phi_prop^2) , log = TRUE)-
             dnorm(h0, mu, sigma/sqrt(1-phi^2) , log = TRUE)

  # the priors of mu and phi
  ar_prob <- ar_prob +
  dnorm(gamma_prop, prior_spec$mu_normal_mean*(1-phi_prop), prior_spec$mu_normal_sd*(1-phi_prop) , log =TRUE) - #prior for gamma
  dnorm(mu*(1-phi), prior_spec$mu_normal_mean*(1-phi), prior_spec$mu_normal_sd*(1-phi) , log =TRUE) +
  dbeta((phi_prop+1)/2, prior_spec$phi_beta_alpha, prior_spec$phi_beta_beta , log =TRUE) -     # transformed phi has a beta
  dbeta((phi+1)/2, prior_spec$phi_beta_alpha, prior_spec$phi_beta_beta , log = TRUE)

  # the proposal densities
  ar_prob <- ar_prob +
  dnorm(phi, 0, sigma * proposal_phi_sd , log =TRUE) -   #phi
  dnorm(phi_prop, 0, sigma_prop * proposal_phi_sd , log =TRUE) +
  dnorm(mu*(1-phi), 0, sigma * proposal_gamma_sd) -               # old gamma proposal
  dnorm(gamma_prop, 0, sigma_prop * proposal_gamma_sd)            # new gamma proposal

  # Now sigma:
  # Proposal is InvGamma(-0.5, 0)
  # Target is Gamma(.5, 1/(2*Bsigma))
  if(shrinksv == "volvol"){
  ar_prob <- ar_prob +
  #(sigma^2 - sigma_prop^2) / (2 * 0.5/ prior_spec$sigma2_gamma_rate) #old
  (sigma^2 - sigma_prop^2)*prior_spec$sigma2_gamma_rate
  #(sigma^2 - sigma_prop^2) / (2 * prior_spec$sigma2_gamma_rate)
  
  # shirnking stationary volatility:
  }else if(shrinksv == "statvolvol"){
  # simplified version:
  # ar_prob <- ar_prob +
  #   0.5*(log(1-phi^2)-log(1-phi_prop^2))-
  #   sigma_prop^2/((1-phi_prop^2)*prior_spec$sigma2_gamma_rate)+
  #   sigma^2/((1-phi^2)*prior_spec$sigma2_gamma_rate)

    
    
    # ar_prob <- ar_prob +
    #   log(sigma_prop)-log(sigma)+
    #   dgamma(sigma_prop^2, 1/2 ,1/((1-phi_prop^2)*prior_spec$sigma2_gamma_rate), log =TRUE)-
    #   dgamma(sigma^2, 1/2 , 1/((1-phi^2)*prior_spec$sigma2_gamma_rate), log =TRUE)
    
    # 02.08.2022 wrong- 2 missing in the denominator, corrected:
    
    ar_prob <- ar_prob +
    log(sigma_prop)-log(sigma)+
    dgamma(sigma_prop^2, 1/2 ,1/(2*(1-phi_prop^2)*prior_spec$sigma2_gamma_rate), log =TRUE)-
    dgamma(sigma^2, 1/2 , 1/(2*(1-phi^2)*prior_spec$sigma2_gamma_rate), log =TRUE)

    
    
    
  }

  # accept/reject
  if(log(runif(1))  < ar_prob){
    return(list(mu = mu_prop, phi = phi_prop, sigma = sigma_prop))
  } else {
    return(list(mu = mu, phi = phi, sigma = sigma))
  }

  }else{ #2block sampler:
    
    shrinksv <- prior_spec$shrinksv
    # calculate bT and BT:
    X <- matrix(c(rep(1,Ttot),c(h0,head(h,-1) )), Ttot , 2)
    gamma_varinv <- 1/prior_spec$gamma_var
    phi_varinv <- 1/prior_spec$phi_var
    B0_inv <- diag(c(gamma_varinv,phi_varinv))
    # BT_inv_11 <- Ttot-1 + B0_inv[1,1]
    # BT_inv_12 <- sum(X[,2])
    # BT_inv_22 <- sum(X[,2]*X[,2])
    # BT <- 1/(BT_inv_11*BT_inv_22 - BT_inv_12^2) * matrix(c(BT_inv_22,-BT_inv_12,-BT_inv_12,BT_inv_11), 2,2, byrow = TRUE)
    
    # old:
    BT <- solve(crossprod(X)+B0_inv)
    bT <- BT%*%t(X)%*%h
    
    gamma <- mu*(1-phi)
    
    # propose sigma using the GIG distribution
    lambda_sig2 <- -Ttot/2
    chi_sig2 <- (1-phi^2)*(h0-mu)^2 + sum( ((h-mu) - phi*(c(h0,head(h,-1))-mu))^2)
    #chi_sig2 <- (1-phi^2)*(h0-mu)^2 + sum( ((h-gamma) - phi*(c(h0,head(h,-1))))^2)
    # chi_sig2 <- (1-phi^2)*(h0-mu)^2+((h[1]-mu)-phi*(h0-mu))^2
    # for(t in 2:Ttot){
    #   chi_sig2 <- chi_sig2 + ((h[t]-mu)-phi*(h[t-1]-mu))^2
    # }
    psi_sig2 <- 2*prior_spec$sigma2_gamma_rate
    sigma_prop <- GIGrvg::rgig(n=1, lambda=lambda_sig2, chi=chi_sig2, psi=psi_sig2)
    
    # not working yet:
    # try using the old method
    cT <- (Ttot)/2
    CT <- .5*((1-phi^2)*(h0-mu)^2 + sum( ((h-gamma) - phi*(c(h0,head(h,-1))))^2))
    sigma_prop <- sqrt(1/rgamma(1,cT, CT))
    
    # Calculate Cholesky decomposition of sig2BT- covariance matrix of gamma and phi:
    chol11 <- sqrt(BT[1,1])*sigma
    chol12 <- (BT[1,2]/chol11)*sigma
    chol22 <- sqrt(BT[2,2]-chol12^2)*sigma
    
    # propose phi using inverse transform sampling from a truncated normal distribution:
    quant_low <-  pnorm(-1, bT[2], chol22, lower.tail = TRUE ,log.p = FALSE)  # use chol11 = sqrt(BT[1,1]*sigma^2)
    quant_high <- pnorm(1, bT[2], chol22,  lower.tail = TRUE ,log.p = FALSE)
    phi_prop <- qnorm((quant_low + runif(1)*(quant_high-quant_low)), bT[2],sd= chol22,lower.tail = TRUE ,log.p = FALSE)
    
    # propose gamma
    gamma_prop <- rnorm(1, bT[1] + chol12*((phi_prop-bT[2])/chol22), sd = chol11) # rewrite in terms of BT?
    
    # Calculate acceptance rates
    # Start with gamma and phi
    mu_prop <- gamma_prop/(1-phi_prop)
    #phi_const <- 1-phi
    proposal_gamma_sd <- sqrt(prior_spec$gamma_var)
    proposal_phi_sd <- sqrt(prior_spec$phi_var)
    
    # first the likelihood
    ar_prob <- dnorm(h0, mu_prop, sigma_prop/sqrt(1-phi_prop^2) , log = TRUE)-
      dnorm(h0, mu, sigma/sqrt(1-phi^2) , log = TRUE)
    
    # the priors of mu and phi
    ar_prob <- ar_prob +
      dnorm(gamma_prop, prior_spec$mu_normal_mean*(1-phi_prop), prior_spec$mu_normal_sd*(1-phi_prop) , log =TRUE) - #prior for gamma
      dnorm(mu*(1-phi), prior_spec$mu_normal_mean*(1-phi), prior_spec$mu_normal_sd*(1-phi) , log =TRUE) +
      dbeta((phi_prop+1)/2, prior_spec$phi_beta_alpha, prior_spec$phi_beta_beta , log =TRUE) -     # transformed phi has a beta
      dbeta((phi+1)/2, prior_spec$phi_beta_alpha, prior_spec$phi_beta_beta , log = TRUE)
    
    # the proposal densities
    ar_prob <- ar_prob +
      dnorm(phi, 0, sigma * proposal_phi_sd , log =TRUE) -   #phi
      dnorm(phi_prop, 0, sigma_prop * proposal_phi_sd , log =TRUE) +
      dnorm(mu*(1-phi), 0, sigma * proposal_gamma_sd) -               # old gamma proposal
      dnorm(gamma_prop, 0, sigma_prop * proposal_gamma_sd)            # new gamma proposal
    
    # Now sigma:
    # Proposal is InvGamma(-0.5, 0)
    # Target is Gamma(.5, 1/(2*Bsigma))
    if(shrinksv == "volvol"){
      ar_prob_sigma <-  (sigma^2 - sigma_prop^2)*prior_spec$sigma2_gamma_rate
     
      # shrinking stationary volatility:
    }else if(shrinksv == "statvolvol"){
      # simplified version:
      # ar_prob <- ar_prob +
      #   0.5*(log(1-phi^2)-log(1-phi_prop^2))-
      #   sigma_prop^2/((1-phi_prop^2)*prior_spec$sigma2_gamma_rate)+
      #   sigma^2/((1-phi^2)*prior_spec$sigma2_gamma_rate)
      
      # ar_prob <- ar_prob +
      #   log(sigma_prop)-log(sigma)+
      #   dgamma(sigma_prop^2, 1/2 ,1/((1-phi_prop^2)*prior_spec$sigma2_gamma_rate), log =TRUE)-
      #   dgamma(sigma^2, 1/2 , 1/((1-phi^2)*prior_spec$sigma2_gamma_rate), log =TRUE)
      
      # 02.08.2022 wrong- 2 missing in the denominator, corrected:
      
      ar_prob_sigma <- 
        log(sigma_prop)-log(sigma)+
        dgamma(sigma_prop^2, 1/2 ,1/(2*(1-phi_prop^2)*prior_spec$sigma2_gamma_rate), log =TRUE)-
        dgamma(sigma^2, 1/2 , 1/(2*(1-phi^2)*prior_spec$sigma2_gamma_rate), log =TRUE)
      
    }
    
    
    # accept/reject
    
    if(log(runif(1))  < ar_prob_sigma){
      sigma <- sigma_prop
    } else {
      # keep sigma to the old value, return always sigma
    }
    
    if(log(runif(1))  < ar_prob){
      return(list(mu = mu_prop, phi = phi_prop, sigma = sigma))
    } else {
      return(list(mu = mu, phi = phi, sigma = sigma))
    }

  } # end of 2-Block sampler

} # end of function



