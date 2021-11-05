library(tidyverse)
library(ggmcmc)
library(nimble)
library(coda)
library(here)


load(here('data','Phoques.Rdata'))
mythin=1
myburn=500
mynit=myburn+mythin*1000


# augmentation des data avec 200 phoque fictif
yAug <- rbind(obs_phoque,matrix(0,ncol = ncol(obs_phoque),nrow = 200))
wtAug <- c(wt_phoque,rep(NA,200))

# prep data  ---------------------
data<-list(yaug=as.matrix(yAug),wt=wtAug)
const=list(nt = ncol(yAug),
           M=nrow(yAug))


# Nimble model   ---------
# M0  -------
M0 <- nimbleCode({
    # Priors
    omega ~ dunif(0, 1)
    p ~ dunif(0, 1)

    # Likelihood
    for (i in 1:M){
        z[i] ~ dbern(omega)

        for (t in 1:nt){
            yaug[i,t] ~ dbern(p.eff[i,t])
            p.eff[i,t] <- z[i] * p
        } #j
    } #i

    # Derived quantities
    N <- sum(z[1:M])
})


# Mt ----------------------------------------------------------------------


Mt1 <- nimbleCode({
    # Priors
    omega ~ dunif(0, 1)
    for(t in 1:nt){
        #fix
        p[t] ~ dunif(0, 1)
    }

    # Likelihood
    for (i in 1:M){
        z[i] ~ dbern(omega)
        for (t in 1:nt){
            yaug[i,t] ~ dbern(p.eff[i,t])
            p.eff[i,t] <- z[i] * p[t]
        } #j
    } #i

    # Derived quantities
    N <- sum(z[1:M])
})
Mt2 <- nimbleCode({
    # Priors
    omega ~ dunif(0, 1)
    mu.p~dlogis(0,1)
    sd.p~dunif(0,5)
    for(t in 1:nt){
        #random
        logit(p[t]) ~ dnorm(mu.p, sd=sd.p)
    }

    # Likelihood
    for (i in 1:M){
        z[i] ~ dbern(omega)
        for (t in 1:nt){
            yaug[i,t] ~ dbern(p.eff[i,t])
            p.eff[i,t] <- z[i] * p[t]
        } #j
    } #i

    # Derived quantities
    N <- sum(z[1:M])
})




# Mb ----------------------------------------------------------------------


Mb <- nimbleCode({
    # Priors
    omega ~ dunif(0, 1)
    p ~ dunif(0, 1)
    c ~ dunif(0, 1)
    # Likelihood
    for (i in 1:M){
        z[i] ~ dbern(omega)

        # First occasion
        yaug[i,1] ~ dbern(p.eff[i,1])
        p.eff[i,1] <- z[i] * p

        # All subsequent occasions
        for (t in 2:nt){
            yaug[i,t] ~ dbern(p.eff[i,t])
            p.eff[i,t] <- z[i] * ( (1-yaug[i,(t-1)]) * p + yaug[i,(t-1)] * c )
        } #j
    } #i

    # Derived quantities
    N <- sum(z[1:M])
})

# Mh ----------------------------------------------------------------------

Mth_x <- nimbleCode({
    # Priors
    omega ~ dunif(0, 1)
    mu.p ~ dlogis(0, 1)
    sd.id~dunif(0,5)
    sd.t~dunif(0,5)
    b_wt~dnorm(0,0.01)
    mu.wt~dnorm(0,0.01)
    sd.wt~dunif(0,3)

    for(t in 1:nt){
        ranef_t[t]~dnorm(0,sd=sd.t)
    }

    for(i in 1:M){
        ranef_id[i]~dnorm(0,sd=sd.id)
        wt[i]~T(dnorm(mu.wt,sd=sd.wt),-6,6)
        for(t in 1:nt){
            logit(p[i,t]) <- mu.p+
                wt[i]*b_wt+
                ranef_id[i]+
                ranef_t[t]
        }
    }


    # Likelihood
    for (i in 1:M){
        z[i] ~ dbern(omega)
        # All subsequent occasions
        for (t in 1:nt){
            yaug[i,t] ~ dbern(p.eff[i,t])
            p.eff[i,t] <- z[i] * p[i,t]
        } #j
    } #i

    # Derived quantities
    N <- sum(z[1:M])
})




# Initial values
inits <- function() list(z = rep(1,nrow(yAug)),
                         b_wt=rnorm(1,0,0.1)
                         )


plot(samples$samples[, "N"])
samples$summary$all.chains %>% head()
samples <- nimbleMCMC(
    code = .x,
    data=data,constants =const,
    inits = inits,
    monitors = c("N",  "omega"#, "z"
                 # ,'b_wt',"mu.p",'mu.wt','sd.wt','sd.id','sd.t'
    ),
    nburnin = myburn,thin = mythin, niter = mynit, nchains = 3,
    WAIC=F,summary = T,samplesAsCodaMCMC = T
)

# run all models in a loop
allMod <- list(m0=M0,mt=Mt1,mtr=Mt2,mb=Mb,mx=Mth_x)
allOut <- map(allMod, function(.x){
    samples <- nimbleMCMC(
        code = .x,
        data=data,constants =const,
        inits = inits,
        monitors = c("N",  "omega"#, "z"
                     # ,'b_wt',"mu.p",'mu.wt','sd.wt','sd.id','sd.t'
                     ),
        # niter = 2000, nburnin = 1000, thin = 1,nchains = 1,
        nburnin = myburn,thin = mythin, niter = mynit, nchains = 3,
        WAIC=F,summary = T,samplesAsCodaMCMC = T
        )
    return(samples)
})




# Jolly-Seber -------------------------------------------------------------

JS <- nimbleCode({
    # Priors and constraints
    for (i in 1:M){
        for (t in 1:(nt-1)){
            phi[i,t] <- mean.phi
        } #t
        for (t in 1:nt){
            p[i,t] <- mean.p
        } #t
    } #i

    mean.phi ~ dunif(0, 1)         # Prior for mean survival
    mean.p ~ dunif(0, 1)           # Prior for mean capture
    psi ~ dunif(0, 1)              # Prior for inclusion probability

    # Dirichlet prior for entry probabilities
    for (t in 1:nt){
        beta[t] ~ dgamma(1, 1)
        b[t] <- beta[t] / sum(beta[1:nt])
    }

    # Convert entry probs to conditional entry probs
    nu[1] <- b[1]
    for (t in 2:nt){
        nu[t] <- b[t] / (1-sum(b[1:(t-1)]))
    } #t

    # Likelihood
    for (i in 1:M){
        # First occasion
        # State process
        w[i] ~ dbern(psi)                  # Draw latent inclusion
        z[i,1] ~ dbern(nu[1])
        # Observation process
        mu1[i] <- z[i,1] * p[i,1] * w[i]
        yaug[i,1] ~ dbern(mu1[i])

        # Subsequent occasions
        for (t in 2:nt){
            # State process
            q[i,t-1] <- 1-z[i,t-1]
            mu2[i,t] <- phi[i,t-1] * z[i,t-1] + nu[t] * prod(q[i,1:(t-1)])
            z[i,t] ~ dbern(mu2[i,t])
            # Observation process
            mu3[i,t] <- z[i,t] * p[i,t] * w[i]
            yaug[i,t] ~ dbern(mu3[i,t])
        } #t
    } #i

    # Calculate derived population parameters
    for (i in 1:M){
        for (t in 1:nt){
            u[i,t] <- z[i,t]*w[i]     # Deflated latent state (u)
        }
    }
    for (i in 1:M){
        recruit[i,1] <- u[i,1]
        for (t in 2:nt){
            recruit[i,t] <- (1-u[i,t-1]) * u[i,t]
        } #t
    } #i
    for (t in 1:nt){
        N[t] <- sum(u[1:M,t])        # Actual population size
        B[t] <- sum(recruit[1:M,t])  # Number of entries
    } #t
    for (i in 1:M){
        Nind[i] <- sum(u[i,1:nt])
        Nalive[i] <- 1-equals(Nind[i], 0)
    } #i
    Nsuper <- sum(Nalive[1:M])         # Superpopulation size
})

inits <- function() list(
    psi=runif(1,0.2,0.7),
    mean.p=runif(1,0.2,.5),
    mean.phi=runif(1,0.2,.5),
    z = matrix(1,nrow(yAug),ncol(yAug)),
                         w=rep(1,nrow(yAug))
                         )

start <- Sys.time()
JS_out <- nimbleMCMC(
    code = JS,
    data=data,constants =const,
    inits = inits,
    monitors =c("psi", "mean.p", "mean.phi", "b", "Nsuper", "N", "B", "nu"),
    # nburnin = myburn,thin = mythin, niter = mynit, nchains = 3,
    niter = 5000, nburnin = 3000, thin = 2,nchains = 2,
    WAIC=F,summary = T,samplesAsCodaMCMC = T
)
dur= Sys.time()-start

plot(JS_out$samples[,'Nsuper'])
