library(rstan) 
library(TruncatedNormal)
#We build a data set in order to check the accuracy of our stan code 
#================================================================================
# a function which randomly generates X and y given mu_alpha, sigma_gamma,
#mu_beta, sigma_beta, mu_gamma, sigma_gamma, lambda and delta
gen_y <- function(n, p, mu_gamma, sigma_gamma, mu_alpha, sigma_alpha, mu_beta, sigma_beta, lambda, delta){
  gamma <- rnorm(p, mu_gamma, sigma_gamma)
  alpha <- rtmvnorm(1, lb = rep(0, 3), mu = mu_alpha, sigma = sigma_alpha)
  beta <- rtmvnorm(1, lb = 0, mu = mu_beta, sigma = sigma_beta)
  X <- matrix(data = rnorm(n*p, mean = 0, sd = 2), nrow = n, ncol = p)
  nu <- rexp(1, lambda) + delta
  rho <- (nu-2)/nu
  #we build recursively h and u
  h <- c(alpha[1])#we assume that u[0]=0
  u <- c(rt(1, nu)* sqrt(rho*h[1]))
  if (u[1] >= 0){ 
    h <- c(h, alpha[1] + alpha[2]*u[1]^2+beta*h[1])
  }else{ 
    h <- append(h, alpha[1] + alpha[3]*u[1]^2+beta*h[1])
  }
  
  for (t in 2:n){
    u <- append(u,rt(1, nu)* sqrt(rho*h[t]))
    if (u[t] >= 0){
      h <- append(h, alpha[1] + alpha[2]*u[t]^2+beta*h[t])
    }else{
      h <- append(h, alpha[1] + alpha[3]*u[t]^2+beta*h[t])
    }
  }
  #we build y
  y <- as.vector(X %*% gamma + u)
  return(list("y"=y, "X"=X, "alpha"=alpha, "beta"=beta, "gamma"=gamma, "nu"=nu))
}



GJR11 <- function(y, X, mu_alpha=rep(1,3), sigma_alpha=diag(x=rep(0.01,3, nrow=3, ncol=3)),
                  mu_beta=1, sigma_beta=0.01,
                  mu_gamma=rep(1,ncol(X)), sigma_gamma=diag(x=rep(0.01,ncol(X)), nrow=ncol(X), ncol=ncol(X)),
                  lambda=1, delta=3, student=FALSE,
                  warmup = is.integer(iter/2), iter=3000, chains=1, cores=1, init=0.1){
  if(student){
    data = list(n=nrow(X),p=ncol(X),y=y, X=X, 
                mu_alpha=mu_alpha, sigma_alpha=sigma_alpha,
                mu_beta=mu_beta, sigma_beta=sigma_beta,
                mu_gamma=mu_gamma, sigma_gamma=sigma_gamma,
                lambda=lambda, delta=delta) 
    fit = stan(file = "models/student-t-GJR11.stan", data = data, iter=iter, chains=chains,
               cores=cores, init=init, control = list(adapt_delta = 0.99))
  }
  else{
    data = list(n=nrow(X),p=ncol(X),y=y, X=X, 
                mu_alpha=mu_alpha, sigma_alpha=sigma_alpha,
                mu_beta=mu_beta, sigma_beta=sigma_beta,
                mu_gamma=mu_gamma, sigma_gamma=sigma_gamma) 
    fit = stan(file = "models/normal-GJR11.stan", data = data, iter=iter, chains=chains,
               cores=cores, init=init, warmup=warmup,control = list(adapt_delta = 0.99))
  }
  return(fit)
}