# nimble intro
library(tidyverse)
library(nimble)
library(coda)
library(ggmcmc)

data <- list(x = c(8, 1, 6, 2, 3, 6, 5, 7, 6, 9),
             y= c(77, 24, 44, 35, 10, 54, 72, 85, 75, 66)  )
const=list(N=10)

code <- nimbleCode({
  beta[1] ~ dnorm(0, sd = 100)
  beta[2] ~ dnorm(0,sd=100)
  sigma ~ dunif(0, 100)

  for(i in 1:N){
    yhat[i] <- beta[1]+beta[2]*x[i]
    y[i] ~ dnorm(yhat[i], sd = sigma)
  }
})

initsFunction <- function() list(beta = rnorm(2,0,1),sigma = runif(1,0,10))

mcmc.out <- nimbleMCMC(code = code,
                       data = data, constants = const,inits = initsFunction,
                       monitors = c("beta", "sigma"),
                       summary=T,samplesAsCodaMCMC = T,
                       thin = 2, niter = 2000, nchains = 3)

# look at output
g_mcmc <- ggs(mcmc.out$samples)
ggs_traceplot(g_mcmc,family='beta')
mcmc.out$summary$all.chains
effectiveSize(mcmc.out$samples)
gelman.diag(mcmc.out$samples)

# compare prior and posterior
priors <- data.frame(
  Parameter=rep(c('beta[1]','beta[2]'),each=3000),
  value=c(rnorm(3000,0,sd=100),rnorm(3000,0,sd=100)),
  prior='prior'
)
g_mcmc %>% filter(Parameter %in% c('beta[1]','beta[2]')) %>%
  select(Parameter,value) %>% mutate(prior='posterior') %>%
  bind_rows(priors) %>%
  ggplot(aes(x=value,fill=prior))+geom_density(alpha=0.5)+facet_wrap(~Parameter)



# Nimble workflow : version longue ----------------------------------------
# build model
mymodel <- nimbleModel(code = code,
           data = data, constants = const,inits = initsFunction())
# Configure MCMC
mymcmcConf <- configureMCMC(mymodel)
mymcmcConf$addMonitors(c("beta", "sigma"))

myMCMC <- buildMCMC(mymcmcConf)

#compile to C++
Cmod <- compileNimble(mymodel)
CMCMC <- compileNimble(myMCMC, project = mymodel)

# run MCMC and get samples
out <- runMCMC(CMCMC,
               summary=T,samplesAsCodaMCMC = T,
                thin = 2,niter = 2000, nchains = 3)
