#Seting data====================================================================
#The historic of these exchange rates comes from https://www.macrotrends.net
USD_JPY <- read.csv("data/dollar-yen-exchange-rate-historical-chart.csv", skip=13)
EUR_USD <- read.csv("data/euro-dollar-exchange-rate-historical-chart.csv", skip=13)
GBP_USD <- read.csv("data/pound-dollar-exchange-rate-historical-chart.csv", skip=13)
names(USD_JPY) <- c("date", "USD_JPY")
names(EUR_USD) <- c("date", "EUR_USD")
names(GBP_USD) <- c("date", "GBP_USD")

#defining the proper year format
USD_JPY$date <- format(as.Date(USD_JPY$date), "%Y-%m-%d")
EUR_USD$date <- format(as.Date(EUR_USD$date), "%Y-%m-%d")
GBP_USD$date <- format(as.Date(GBP_USD$date), "%Y-%m-%d")

#taking the prices of May and June 2020
USD_JPY <- USD_JPY[USD_JPY$date>=as.Date("2020-06-20") & USD_JPY$date<=as.Date("2020-06-30"),]
EUR_USD <- EUR_USD[EUR_USD$date>=as.Date("2020-06-20") & EUR_USD$date<=as.Date("2020-06-30"),]
GBP_USD <- GBP_USD[GBP_USD$date>=as.Date("2020-06-20") & GBP_USD$date<=as.Date("2020-06-30"),]
#we notice that the markets are closed the first and second days of each January


#We have the rate USD/JPY, we want JPY/USD
inver <- function(x){return(1/x)}
JPY_USD <- USD_JPY
JPY_USD$USD_JPY <- apply(USD_JPY[2],2, FUN=inver)
names(JPY_USD) <- c("date", "JPY_USD")

#We merge all the data in one data frame for convenience
rates <- merge(EUR_USD, GBP_USD, by = "date")
rates <- merge(rates, JPY_USD, by = "date")
names(rates) <- c("date", "EUR_USD", "GBP_USD", "JPY_USD")

rm("EUR_USD", "GBP_USD", "JPY_USD", "USD_JPY", "inver")

#Stating relevant hyperparameters===============================================
#we carry out a linear regression with bootstrap to evaluate the mean and 
#covariance matrix of gamma
p <- 20
mu_gamma_b <- c()
for(i in 1:p){
  data <- rates[sample(x=nrow(rates)
                         ,replace=TRUE),c("EUR_USD", "GBP_USD", "JPY_USD")]
  linearMod <- lm(EUR_USD ~ GBP_USD + JPY_USD, data=data)
  mu_gamma_b <- rbind(mu_gamma_b, linearMod$coefficients)
}
mu_gamma <- c(mean(mu_gamma_b[,2]), mean(mu_gamma_b[,3]))#we don’t want "intercept which the constant

sigma_gamma <- var(mu_gamma_b)[2:nrow(var(mu_gamma_b)),2:ncol(var(mu_gamma_b))]
sigma_gamma <- diag(diag(sigma_gamma), nrow=2, ncol=2)
sigma_gamma <- sqrt(sigma_gamma)
#we neglect the covariace between two diffirents coefs

#for alpha and beta, we fix mu=X*gamma_hat. Then we evaluate the most likely 
#value of alpha_hat and beta_hat by maximizing the likelihood of a GJR model 
#with bootstrap. To do so, we use the package "rugarch"
library(rugarch)

spec <-  ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
                  variance.model = list(model = "gjrGARCH"), distribution.model = "sstd")
fit <-  ugarchfit(spec, data = rates$EUR_USD)
#by watching the results of the bootstrap regression, we state 
mu_beta <- 0.9
sigma_beta <- 0.05
mu_alpha <- as.vector(c(1,1,1))
sigma_alpha <- diag(x=rep(0.04,3),nrow=3,ncol=3)

#Running stan model=============================================================
source("GJR11.R") 
y <- as.vector(rates$EUR_USD)
X <- cbind(rates$GBP_USD, rates$JPY_USD)
X <- cbind(X)#we add a column of 1 to enable the existance of a constant
#normal GJR
beg <- Sys.time()
fit_normal <- GJR11(y,X, mu_alpha=mu_alpha, sigma_alpha=sigma_alpha, 
             mu_beta=mu_beta, sigma_beta=sigma_beta, 
             mu_gamma=mu_gamma, sigma_gamma=sigma_gamma,
             lambda=2, delta=3, student=FALSE, iter=100000, chains=4, cores=4)
cat("time:", Sys.time()-beg)
summary(fit_normal, c("alpha", "beta", "gamma"))$summary
#student GJR
beg <- Sys.time()
fit_student <- GJR11(y,X, mu_alpha=mu_alpha, sigma_alpha=sigma_alpha, 
                    mu_beta=mu_beta, sigma_beta=sigma_beta, 
                    mu_gamma=mu_gamma, sigma_gamma=sigma_gamma,
                    lambda=2, delta=3, student=TRUE, iter=100000, chains=4, cores=4,
                    warmup = as.integer(iter*(3/4)))
cat("time:", Sys.time()-beg)
summary(fit_student, c("alpha", "beta", "gamma", "nu"))$summary

good_model <- fit_student
round( summary(good_model, c("alpha", "beta", "gamma", "nu"))$summary  , 4)
#autocorrelation
acf(x = rates[2:4,2:4], lag.max = 2,main=NA)
#some plots
library(bayesplot)
posterior <- as.matrix(good_model)
mcmc_areas(posterior, pars = c("alpha"))

library("dplyr")
color_scheme_set("brightblue")
good_model %>%
  posterior_predict(draws = 500) %>%
  ppc_stat_grouped(y = mtcars$mpg,
                   group = mtcars$carb,
                   stat = "median")

#Bayes factor, comparison with a centered normal priors
library(BayesFactor)
extractBF(good_model)

library(brms)
library(bridgesampling)
bf(optimized_hyperparam, normal_hyperparam)



library(bridgesampling)
optimized_hyperparam <- good_model
normal_hyperparam <- GJR11(y,X, student=TRUE, iter=100000, chains=4, cores=4,
                           warmup = as.integer(iter*(3/4)))
summary(normal_hyperparam, c("alpha", "beta", "gamma", "nu"))$summary

optimized_hyperparam_2 <- bridge_sampler(optimized_hyperparam, silent = TRUE)
normal_hyperparam_2 <- bridge_sampler(normal_hyperparam, silent = TRUE)
BF_att <- bridgesampling::bf(optimized_hyperparam_2, normal_hyperparam_2)

BF_att <- bridgesampling::bf(optimized_hyperparam, normal_hyperparam)

#================================================================================
#Now we take into account only the exchange rate GBP/USD
p <- 20
mu_gamma_b <- c()
for(i in 1:p){
  data <- rates[sample(x=nrow(rates)
                       ,replace=TRUE),c("EUR_USD", "GBP_USD")]
  linearMod <- lm(EUR_USD ~ GBP_USD, data=data)
  mu_gamma_b <- rbind(mu_gamma_b, linearMod$coefficients)
}
mu_gamma <- c(mean(mu_gamma_b[,1]), mean(mu_gamma_b[,2]))

sigma_gamma <- var(mu_gamma_b)
sigma_gamma <- diag(diag(sigma_gamma), nrow=2, ncol=2)
sigma_gamma <- sqrt(sigma_gamma)
#we neglect the covariace between two diffirents coefs

#for alpha and beta, we fix mu=X*gamma_hat. Then we evaluate the most likely 
#value of alpha_hat and beta_hat by maximizing the likelihood of a GJR model 
#with bootstrap. To do so, we use the package "rugarch"
library(rugarch)

spec <-  ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
                    variance.model = list(model = "gjrGARCH"), distribution.model = "sstd")
fit <-  ugarchfit(spec, data = rates$EUR_USD)
#by watching the results of the bootstrap regression, we state 
mu_beta <- 0.9
sigma_beta <- 0.05
mu_alpha <- as.vector(c(1,1,1))
sigma_alpha <- diag(x=rep(0.04,3),nrow=3,ncol=3)
#
iter <- 1000
y <- as.vector(rates$EUR_USD)
X <- cbind(rep(1, length(rates$GBP_USD)), rates$GBP_USD)#we add a column of 1 to enable the existance of a constant
beg <- Sys.time()
fit_normal_2 <- GJR11(y,X, mu_alpha=mu_alpha, sigma_alpha=sigma_alpha, 
                    mu_beta=mu_beta, sigma_beta=sigma_beta, 
                    mu_gamma=mu_gamma, sigma_gamma=sigma_gamma,
                    student=FALSE, iter=iter, chains=4, cores=4, 
                    warmup = as.integer(iter*(3/4)), init=1)
cat("time:", Sys.time()-beg)
round(summary(fit_normal_2, c("alpha", "beta", "gamma"))$summary, 3)

beg <- Sys.time()
fit_student_2 <- GJR11(y,X, mu_alpha=mu_alpha, sigma_alpha=sigma_alpha, 
                    mu_beta=mu_beta, sigma_beta=sigma_beta, 
                    mu_gamma=mu_gamma, sigma_gamma=sigma_gamma,
                    student=TRUE, iter=iter, chains=4, cores=4, 
                    warmup = as.integer(iter*(3/4)), init=1)
cat("time:", Sys.time()-beg)
round(summary(fit_normal_2, c("alpha", "beta", "gamma", "nu"))$summary, 3)
