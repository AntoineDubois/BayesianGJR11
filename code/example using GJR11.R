source("GJR11.R")

n <- 50
p <- 4
mu_gamma <- as.vector(rep(1,p))
sigma_gamma <- diag(x = rep(0.01,p), nrow = p, ncol = p)
mu_alpha <- as.vector(rep(0.1,3))
sigma_alpha <- diag(x = rep(0.05,3), nrow = 3, ncol = 3)
mu_beta <- 1
sigma_beta <- 0.05

lambda <- 2
delta <- 3

A <- gen_y(n, p, mu_gamma, sigma_gamma, mu_alpha, sigma_alpha, mu_beta, sigma_beta, lambda, delta)
y <- A$y
X <- A$X

#normal GJR
fit <- GJR11(y,X)
summary(fit, c("alpha", "beta", "gamma"))$summary
cat("alpha:",A$alpha,"beta:",A$beta,"gamma:",A$gamma,"nu:",A$nu)
#student GJR
fit <- GJR11(y,X, student = TRUE)
summary(fit, c("alpha", "beta", "gamma", "nu"))$summary
cat("alpha:",A$alpha,"beta:",A$beta,"gamma:",A$gamma,"nu:",A$nu)