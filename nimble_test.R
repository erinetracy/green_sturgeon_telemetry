#Nimble test

#i installed Rtools 4.5 previous to this
install.packages("nimble")
library(nimble)

#test if nimble is working
code <- nimbleCode({ x ~ dnorm(0, 1) })
model <- nimbleModel(code)
cModel <- compileNimble(model)

#install more packages
install.packages(c("nimbleHMC", "mcmcplots", "nimbleEcology", "compareMCMCs"))
# says mcmcplots isnt available for this version of R, everything else is fine
library(nimbleEcology)
library(nimbleHMC)
library(compareMCMCs)
install.packages(c("CARBayesdata","sp","spdep","classInt", "glmmTMB"))
library(CARBayesdata)
library(sp)
library(spdep)
library(classInt)
library(glmmTMB)


#lesson 0
data(cars)
fit <- lm(dist ~ speed, data = cars)
summary(fit)

library(nimble)

modelcode <- nimbleCode({
  ## Priors
  beta[1] ~ dnorm(0,10) 
  beta[2] ~ dnorm(0,10)
  sigma ~ dunif(0,20)
  
  ## Build the Likelihood
  for (i in 1:n){
    mu[i] <- beta[1] + beta[2]*speed[i]
    dist[i] ~ dnorm(mean = mu[i], sd = sigma)
  }
})

model_constants <- list(speed = cars$speed, n = nrow(cars))
model_data <- list(dist = cars$dist)

model <- nimbleModel(code = modelcode, constants = model_constants,
                     data = model_data)
## model$getNodeNames(determOnly = TRUE)
## model$getNodeNames(stochOnly = TRUE) ## Data and Priors
## model$getNodeNames(stochOnly = TRUE, includeData = FALSE) ## Priors

## Based on lm fit:
model$beta <- as.numeric(coef(fit))
## Adjust MLE for bias correction to match things below.
model$sigma <- sigma(fit)*(nrow(cars)-1)/nrow(cars)

## Calculate the log posterior density (and update all the parameters)
model$calculate()
## Just get the log likelihood.
model$getLogProb(model$getNodeNames(dataOnly = TRUE))
## Identical to what was fit by lm.
logLik(fit)

cmodel <- compileNimble(model)
cmodel$calculate()
microbenchmark::microbenchmark(
  Rversion = model$calculate(),
  Cversion = cmodel$calculate()
)

mcmc.out <- nimbleMCMC(code = modelcode, constants = model_constants,
                       data = model_data,
                       nchains = 3, niter = 10000, 
                       inits = \(){list(beta = rnorm(2), sigma = runif(1, 1,4))},
                       nburnin = 3000, samplesAsCodaMCMC = TRUE,
                       summary = TRUE, monitors = c('beta','sigma'))

mcmc.out$summary


