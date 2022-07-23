draw_h <- function(data, mu, phi, sigma, r,mix_mean, mix_varinv){


  # Define variables
  sigma2 <- sigma^2
  sigma2inv = 1/sigma2
  #Matrix form
  #mix_varinv_diag <- as(diag(mix_varinv), "sparseMatrix")    #bandSparse( Ttot, k = 0, diag = list(mix_varinv))
  #sigma2inv_diag <-  as(diag(rep(sigma2inv, Ttot)), "sparseMatrix")     #bandSparse( Ttot, k = 0, diag = list(sigma2inv))
  # vector form
  #mix_varinv_diag <- mix_varinv
  #sigma2inv_diag <-  rep(sigma2inv, Ttot)


  # Create the covariance matrix omega

  Omega_diag <- mix_varinv[r[-1]] +(1+phi^2)*sigma2inv  # taken from stochvol
  #Omega_diag <- mix_varinv_diag[r] +(1+phi^2)*sigma2inv_diag
  Omega_offdiag <- -phi*sigma2inv             #scalar for now
  Omega_diag[1] =  sigma2inv    # taken from stochvol
  #Omega_diag[1] <- mix_varinv[r[1]] +sigma2inv
  Omega_diag[Ttot] <- mix_varinv[r[Ttot-1]] +sigma2inv   # taken from stochvol
  # Create the covector c
  covector <- (data[-1] - mix_mean[r[-1]])*mix_varinv[r[-1]] + mu*(1-phi)^2*sigma2inv
  covector[1] = mu * (1 - phi) * sigma2inv # taken from stochvol
  #covector[1] <- (data[1] - mix_mean[r[1]])*mix_varinv[r[1]] + mu*(1-phi)*sigma2inv
  covector[Ttot] <- (data[Ttot-1] - mix_mean[r[Ttot-1]])*mix_varinv[r[Ttot-1]] + mu*(1-phi)*sigma2inv

  # double check r-index

  # Compute the Cholesky decomposition

  #Matrix::Cholesky(forceSymmetric(Sigma.inv), perm = FALSE, LDL = FALSE)
  #steptwo <- Matrix::solve(L.chol, Omega.offdiag, system = "L")

  #cholesky_tridiagonal
  #old:
  #Matrix form
  # chol_diag <- diag(Ttot)
  # chol_offdiag <- diag(Ttot)
  # # vector form:
  # chol_diag  <- rep(NA, Ttot)
  # chol_offdiag <- rep(NA, Ttot-1)
  # chol_diag[1] <- sqrt(Omega_diag[1])
  # for(i in 2:Ttot){
  # chol_offdiag[i-1] <- Omega_offdiag/chol_diag[i-1]
  # chol_diag[i] <- sqrt(Omega_diag[i]-chol_offdiag[i-1]^2)
  # }

  #new:
  Omega_chol <- cholesky_tridiagonal(Omega_diag, Omega_offdiag)
  
  #old:
  # Sample h using forward and backward back-substitution
  #forward_algorithm
  # htmp <- rep(NA, Ttot-1)
  # htmp[1] <- covector[1]/chol_diag[1]
  # for(i in 2:Ttot){
  #   htmp[i] <- (covector[i] - chol_offdiag[i-1]*htmp[i-1])/chol_diag[i]
  # }
  # 
  # e <- rnorm(Ttot)
  # htmp <- htmp + e
  # 
  # #backward_algorithm
  # h <- rep(NA, Ttot)
  # h[Ttot] <- htmp[Ttot]/chol_diag[Ttot]
  # for(i in (Ttot-1):1){
  #   h[i] <- (htmp[i] - chol_offdiag[i] * h[i+1])/chol_diag[i]
  # }

  #new:
  htmp <- forward_algorithm(Omega_chol$chol_diag, Omega_chol$chol_offdiag, covector)
  e <- rnorm(Ttot)
  htmp <- htmp + e
  h <- backward_algorithm(Omega_chol$chol_diag, Omega_chol$chol_offdiag, htmp)

  return(h)

}
