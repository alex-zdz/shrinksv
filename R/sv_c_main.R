sv_c_main <- function(y, draws, burnin , thinpara ,  prior_spec, startpara, startlatent, expert
                      ){
  
    #set.seed(1)
    #log_data2_normal, mu, phi, sigma, h0, h, r, prior_spec, expert
    Ttot <- length(y)
    data <- log(y^2)
    
    # number of MCMC draws
    niter = burnin + draws

    # initialize the variables
    mu <- startpara$mu
    phi <- startpara$phi
    sigma <- startpara$sigma
    h0 <- startlatent$latent0 
    h <- startlatent$latent      # contains h1 to hT, but not h0!
    r <- startpara$r  # mixture indicators, need to be between 1 and 10 inclusive

    # storage
    para_draws <- draws / thinpara
    para_store <- matrix(NA, para_draws , 3)
    latent0_store <- rep(NA , draws)
    latent_store <- matrix(NA , draws , Ttot)
    r_store <- matrix( NA , draws , Ttot)

    # Constants and their transforms from Omori et al (2007).
    mix_prob <- c(.00609, .04775, .13057, .20674, .22715, .18842, .12047, .05591, .01575, .00115)
    mix_mean <- c(1.92677, 1.34744, .73504, .02266, -.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
    mix_var <- c(.11265, .17788, .26768, .40611, .62699, .98583, 1.57469, 2.54498, 4.16591, 7.33342)
    mix_sd <- c(sqrt(mix_var))
    mix_varinv <- c(1 / mix_var)
    mix_2varinv <- c(0.5 * mix_varinv)
    mix_pre <- log(mix_prob/mix_sd)

    #  progress bar
    pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                         max = niter, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")

    # main mcmc loop
    for( i in 1:niter){
    
      # Step (c): sample indicators
      r <- sample_r(data , h , mix_prob , mix_mean , mix_var)

      # Step (a): sample the latent volatilities h:
      latent_new <- draw_h(data, mu, phi, sigma, r, mix_mean , mix_varinv)
      h0 <- latent_new[1]
      h <- latent_new[-1]
      
      # Omega_ret <- draw_latent_1(data, mu, phi, sigma, r, mix_mean , mix_varinv)
      # Omega_diag <- Omega_ret$Omega_diag
      # Omega_offdiag <- Omega_ret$Omega_offdiag
      # covector <- Omega_ret$covector
      # Omega_chol <- cholesky_tridiagonal(Omega_diag, Omega_offdiag)
      # htmp <- forward_algorithm(Omega_chol$chol_diag, Omega_chol$chol_offdiag, covector)
      # e <- rnorm(Ttot+1)
      # htmp <- htmp + e
      # latent_new <- backward_algorithm(Omega_chol$chol_diag, Omega_chol$chol_offdiag, htmp)
      # h0 <- latent_new[1]
      # h <- latent_new[-1]
      
      # Step (b): sample mu, phi, sigma
      parameter_draw <- sample_theta(mu, phi, sigma, h0, h, prior_spec,expert = expert )
      mu <-  parameter_draw$mu
      phi <- parameter_draw$phi
      sigma <- parameter_draw$sigma
      #print("done sampling para")
      
      # Storing
      if(i > burnin){
      para_store[i-burnin,1] <- parameter_draw$mu
      para_store[i-burnin,2] <- parameter_draw$phi
      para_store[i-burnin,3] <- parameter_draw$sigma

      latent0_store[i-burnin] <- latent_new[1]
      latent_store[i-burnin,] <- latent_new[-1]
      r_store[i-burnin , ] <- r
      }


      setTxtProgressBar(pb, i)

    }  # END main MCMC loop

    # Progress bar
    # if (verbose) {
    #   if (single_chain) {
    #     progressbar_finish(N);  // finalize progress bar
    #   } else {
    #     chain_print_finish(chain);
    #   }
    # }


return(list(para_store = para_store , latent0_store = latent0_store , latent_store = latent_store , r_store = r_store))



}  # end of function
