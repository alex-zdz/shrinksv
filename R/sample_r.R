sample_r <- function( data , h , mix_prob , mix_mean , mix_var ){

  datanorm <- data-h
  #datanorm <- data - c(h0 , h )

  # Constants and their transforms from Omori et al (2007).
  mix_sd <- c(sqrt(mix_var))
  mix_varinv <- c(1 / mix_var)
  mix_2varinv <- c(0.5 * mix_varinv)

  mixprob_post <- mapply(function(mix_prob , mean , sd ) mix_prob * dnorm(datanorm,mean,sd),mix_prob = mix_prob, mean=mix_mean, sd=sqrt(mix_var))

  # Get CDF:
  r_cdf_unnorm <- apply(t(mixprob_post) , 2 , cumsum)
  r_cdf_unnorm[ 10 , ][r_cdf_unnorm[ 10 , ]==0] <- 1
  r_cdf <- apply(r_cdf_unnorm , 1 , function(x) x/r_cdf_unnorm[ 10 , ])  # divide by 0

  # inverse_transform_sampling
  r <- rep(NA , Ttot)
  usamp <- runif(Ttot)
  r <- apply(r_cdf - usamp,1, function(x) min(which(x>0)))

return(r)

 }
