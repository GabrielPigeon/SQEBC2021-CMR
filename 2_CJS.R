
library(tidyverse)
library(ggmcmc)
library(nimble)
library(coda)
library(here)

load(here('data','reindeer.Rdata'))
mythin=1
myburn=500
mynit=myburn+mythin*1000



const=list(n.occasions=ncol(obs),
           nind=nrow(obs),
           f=yrbirth)
data=list(obs=obs,
          z=obs)

# give guessable values for z (real state, alive/dead)
for(i in 1:nrow(obs)){
    ts=range(which(obs[i,]==1))
    data$z[i,ts[1]:ts[2]] <- 1
    if(ts[2]<ncol(obs)) data$z[i,(ts[2]+1):ncol(obs)] <- NA
}


# CJS: c-c  --------
CJS_cc <- nimbleCode({
phi ~ dunif(0, 1)         # Prior for mean survival
p ~ dunif(0, 1)           # Prior for mean recapture
# Likelihood
for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
        # State process
        mu1[i,t] <- phi * z[i,t-1]
        z[i,t] ~ dbern(mu1[i,t])
        # Observation process
        mu2[i,t] <- p * z[i,t]
        obs[i,t] ~ dbern(mu2[i,t])
    } #t
} #i
})




initsFunction <- function(){
    list(phi = runif(1,0.7,0.9),p = runif(1,.4,.7))
}

CJS_cc.out <- nimbleMCMC(code = CJS_cc,
                       data = data, constants = const,inits = initsFunction,
                       monitors = c("phi", "p"),
                       summary=T,samplesAsCodaMCMC = T,WAIC = F,
                       nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)

CJS_cc.out$summary$all.chains

ggout <- CJS_cc.out$samples %>% ggs()

ggout%>%ggs_traceplot()
ggout %>%ggs_caterpillar()


# CJS: t-t  --------
CJS_tt <- nimbleCode({
    for (t in 1:(n.occasions-1)){
        phi[t] ~ dunif(0, 1)        # Priors for time-spec. survival
        p[t] ~ dunif(0, 1)         # Priors for time-spec. recapture
    }
    # Likelihood
    for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
            # State process
            mu1[i,t] <- phi[t-1] * z[i,t-1]
            z[i,t] ~ dbern(mu1[i,t])
            # Observation process
            mu2[i,t] <- p[t-1] * z[i,t]
            obs[i,t] ~ dbern(mu2[i,t])
        } #t
    } #i
})

initsFunction <- function(){
    list(phi = runif(ncol(obs)-1,0.7,0.9),
         p = runif(ncol(obs)-1,.4,.7)
    )
}

start <- Sys.time()
CJS_tt.out <- nimbleMCMC(code = CJS_tt,
                         data = data, constants = const,inits = initsFunction,
                         monitors = c("phi", "p"),
                         summary=T,samplesAsCodaMCMC = T,WAIC = F,
                         nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)
CJS_tt.dur <- Sys.time()-start


CJS_tt.out$summary$all.chains

ggout <-  CJS_tt.out$samples %>% ggs()

ggout %>% filter(Parameter %in% c('phi[4]', 'p[9]'))%>%
    ggs_traceplot

CJS_tt.out$samples %>% ggs() %>% ggs_caterpillar(family = 'phi',sort=F)



# CJS: t-t random effect --------
CJS_ttr <- nimbleCode({
    for (t in 1:(n.occasions-1)){
        logit(phi[t]) <- phi.mu + ran.phi[t]
        ran.phi[t] ~ dnorm(0,sd=phi.sd)
        logit(p[t]) <- p.mu + ran.p[t]
        ran.p[t] ~ dnorm(p.mu, sd=p.sd)
    }
    phi.mu~dlogis(0,1)
    p.mu~dlogis(0,1)
    phi.sd ~ dunif(0,10)
    p.sd ~ dunif(0,10)

    # Likelihood
    for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
            # State process
            mu1[i,t] <- phi[t-1] * z[i,t-1]
            z[i,t] ~ dbern(mu1[i,t])
            # Observation process
            mu2[i,t] <- p[t-1] * z[i,t]
            obs[i,t] ~ dbern(mu2[i,t])
        } #t
    } #i
})

initsFunction <- function(){
    list(ran.phi = rnorm(ncol(obs)-1,0,0.5),
         ran.p = rnorm(ncol(obs)-1,0,.5),
         p.mu= rnorm(1,0,0.5),
         phi.mu=rnorm(1,0,0.5),
         p.sd=runif(1,0.1,2),
         phi.sd=runif(1,0.1,2)
    )
}

CJS_ttr.out <- nimbleMCMC(code = CJS_ttr,
                         data = data, constants = const,inits = initsFunction,
                         monitors = c("phi", 'phi.mu', 'p.mu','ran.phi', 'ran.p', 'phi.sd', 'p.sd'),
                         summary=T,samplesAsCodaMCMC = T,WAIC = F,
                         nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)


CJS_ttr.out$summary$all.chains[c('phi.mu', 'p.mu', 'phi.sd', 'p.sd'),]


ggout <-  CJS_ttr.out$samples %>% ggs() %>%
           filter(!grepl('ran.',Parameter))


ggout %>% filter(Parameter %in% c('phi.mu', 'p.mu', 'phi.sd', 'p.sd')) %>%
    ggs_traceplot()
ggout %>% filter(Parameter %in% c('phi.mu', 'p.mu', 'phi.sd', 'p.sd')) %>%
    ggs_density()

gphi <- ggout %>%  filter(grepl('phi\\[',Parameter)) %>%
    ggs_caterpillar(family = 'phi',sort=F)+xlim(0.6,1)


# CJS: environmental covariate --------
CJS_tcov <- nimbleCode({
    for (t in 1:(n.occasions-1)){
        logit(phi[t]) <- phi.mu + slope*env[t]+ ran.phi[t]
        ran.phi[t] ~ dnorm(0,sd=phi.sd)
        logit(p[t]) <- p.mu + ran.p[t]
        ran.p[t] ~ dnorm(p.mu, sd=p.sd)
    }
    phi.mu~dnorm(0,0.001)
    p.mu~dlogis(0,1)
    phi.sd ~ dunif(0,3)
    p.sd ~ dunif(0,3)
    slope ~ dnorm(0,0.001)

    # Likelihood
    for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
            # State process
            mu1[i,t] <- phi[t-1] * z[i,t-1]
            z[i,t] ~ dbern(mu1[i,t])
            # Observation process
            mu2[i,t] <- p[t-1] * z[i,t]
            obs[i,t] ~ dbern(mu2[i,t])
        } #t
    } #i
})

initsFunction <- function(){
    list(ran.phi = rnorm(ncol(obs)-1,0,0.5),
         ran.p = rnorm(ncol(obs)-1,0,.5),
         p.mu= rnorm(1,0,0.5),
         phi.mu=rnorm(1,0,0.5),
         p.sd=runif(1,0.1,2),
         phi.sd=runif(1,0.1,2),
         slope=rnorm(1,0,0.3)
    )
}

const$env <- ros

CJS_tcov.out <- nimbleMCMC(code = CJS_tcov,
                          data = data, constants = const,inits = initsFunction,
                          monitors = c('phi.mu', 'p.mu','ran.phi', 'ran.p','slope', 'phi.sd', 'p.sd'),
                          summary=T,samplesAsCodaMCMC = T,WAIC = F,
                          nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)



# age ---------------------------------------------------------------------

data$age <- age
ageC <-c(1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
ageC <-c(1, 2, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5,
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)
data$age <- apply(age,1:2,function(x) ageC[x])

const$nbAge <- max(data$age,na.rm = T)

CJS_Age <- nimbleCode({
    for(i in 1:nind){
        for (t in f[i]:(n.occasions-1)){
            logit(phi[i,t]) <- phi.mu[age[i,t]] + ran.phi[t]
        }
    }
    for (t in 1:(n.occasions-1)){
        logit(p[t]) <- p.mu + ran.p[t]
        ran.p[t] ~ dnorm(p.mu, sd=p.sd)
        ran.phi[t] ~ dnorm(0,sd=phi.sd)
    }
    p.mu~dlogis(0,1)
    phi.sd ~ dunif(0,3)
    p.sd ~ dunif(0,3)

    for(a in 1:nbAge){
        phi.mu[a]~dlogis(0,1)
    }

    # Likelihood
    for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
            # State process
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]
            z[i,t] ~ dbern(mu1[i,t])
            # Observation process
            mu2[i,t] <- p[t-1] * z[i,t]
            obs[i,t] ~ dbern(mu2[i,t])
        } #t
    } #i
})

initsFunction <- function(){
    list(ran.phi = rnorm(ncol(obs)-1,0,0.5),
         ran.p = rnorm(ncol(obs)-1,0,.5),
         p.mu= rnorm(1,0,0.5),
         p.sd=runif(1,0.1,2),
         phi.sd=runif(1,0.1,2),
         phi.mu=rnorm(const$nbAge,0,0.5),
         slope=rnorm(1,0,0.3)
    )
}

CJS_Age.out <- nimbleMCMC(code = CJS_Age,
                              data = data, constants = const,inits = initsFunction,
                              monitors = c('phi.mu', 'p.mu','ran.phi', 'ran.p', 'phi.sd', 'p.sd'),
                              summary=T,samplesAsCodaMCMC = T,WAIC = F,
                          nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)



ggout <-  CJS_Age.out$samples %>% ggs() %>%
    filter(!grepl('ran.',Parameter))

ggout %>% mutate(value=ilogit(value)) %>%
    ggs_caterpillar(family = 'phi.mu',sort = F,horizontal = F)+
    scale_y_discrete(labels=c('0','1','2','3-8','9-11','12+'))+
    theme(axis.text.x = element_text(angle = 0))+
    labs(x='survie',y='age')



# random yr age effects ---------------------------------------------------

CJS_AgeR <- nimbleCode({
    p.mu~dlogis(0,1)
    phi.sd ~ dunif(0,3)
    p.sd ~ dunif(0,3)
    Omega[1:nbAge, 1:nbAge] ~ dwish(R[1:nbAge, 1:nbAge], nbAge+1) # Priors for variance-covariance

    for (t in 1:(n.occasions-1)){
        logit(p[t]) <- p.mu + ran.p[t]
        ran.p[t] ~ dnorm(p.mu, sd=p.sd)
        ran.phi[t] ~ dnorm(0,sd=phi.sd)
        eps.yr[t,1:nbAge]~ dmnorm(zero[1:nbAge], Omega[1:nbAge, 1:nbAge])
        for(a in 1:nbAge){
            ranef.yr[t,a] <-  xi[a]*eps.yr[t,a]
        }
    }

    for(a in 1:nbAge){
        phi.mu[a]~dlogis(0,1)
        xi[a]~dunif(0,10)
    }

    for(i in 1:nind){
        for (t in f[i]:(n.occasions-1)){
            logit(phi[i,t]) <- phi.mu[age[i,t]] + ran.phi[t]
        }
    }
    # Likelihood
    for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
            # State process
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]
            z[i,t] ~ dbern(mu1[i,t])
            # Observation process
            mu2[i,t] <- p[t-1] * z[i,t]
            obs[i,t] ~ dbern(mu2[i,t])
        } #t
    } #i


    # derived corr matrix and actual sigma
    Sigma.raw[1:nbAge, 1:nbAge] <- inverse(Omega[1:nbAge, 1:nbAge])
    for (k in 1:nbAge){
        for (k.prime in 1:nbAge){
            rho[k,k.prime] <- Sigma.raw[k,k.prime]/
                sqrt(Sigma.raw[k,k]*Sigma.raw[k.prime,k.prime]) # Correlations
        }
        Sigma[k] <- abs(xi[k])*sqrt(Sigma.raw[k,k])
    }
})

initsFunction <- function(){
    list(ran.phi = rnorm(ncol(obs)-1,0,0.5),
         ran.p = rnorm(ncol(obs)-1,0,.5),
         p.mu= rnorm(1,0,0.5),
         p.sd=runif(1,0.1,2),
         phi.sd=runif(1,0.1,2),
         phi.mu=rnorm(const$nbAge,0,0.5),
         slope=rnorm(1,0,0.3)
    )
}

CJS_Age.out <- nimbleMCMC(code = CJS_Age,
                          data = data, constants = const,inits = initsFunction,
                          monitors = c('phi.mu', 'p.mu','ran.phi', 'ran.p', 'phi.sd', 'p.sd'),
                          summary=T,samplesAsCodaMCMC = T,WAIC = F,
                          nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)

ggout <-  CJS_Age.out$samples %>% ggs() %>%
    filter(!grepl('ran.',Parameter))

ggout %>% mutate(value=ilogit(value)) %>%
    ggs_caterpillar(family = 'phi.mu',sort = F,horizontal = F)+
    scale_y_discrete(labels=c('0','1','2','3-8','9-11','12+'))+
    theme(axis.text.x = element_text(angle = 0))+
    labs(x='survie',y='age')

# cov:age -----------------------------------------------------------------

CJS_tcovAge <- nimbleCode({
    for(i in 1:nind){
        for (t in f[i]:(n.occasions-1)){
            logit(phi[i,t]) <- phi.mu[age[i,t]] +  slope*env[t] + ran.phi[t]
        }
    }
    for (t in 1:(n.occasions-1)){
            logit(p[t]) <- p.mu + ran.p[t]
            ran.p[t] ~ dnorm(p.mu, sd=p.sd)
            ran.phi[t] ~ dnorm(0,sd=phi.sd)
    }
    p.mu~dlogis(0,1)
    phi.sd ~ dunif(0,3)
    p.sd ~ dunif(0,3)
    slope ~ dnorm(0,0.001)

    for(a in 1:nbAge){
    phi.mu[a]~dnorm(0,0.001)
    }

    # Likelihood
    for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
            # State process
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]
            z[i,t] ~ dbern(mu1[i,t])
            # Observation process
            mu2[i,t] <- p[t-1] * z[i,t]
            obs[i,t] ~ dbern(mu2[i,t])
        } #t
    } #i
})

initsFunction <- function(){
    list(ran.phi = rnorm(ncol(obs)-1,0,0.5),
         ran.p = rnorm(ncol(obs)-1,0,.5),
         p.mu= rnorm(1,0,0.5),
         p.sd=runif(1,0.1,2),
         phi.sd=runif(1,0.1,2),
         phi.mu=rnorm(const$nbAge,0,0.5),
         slope=rnorm(1,0,0.3)
    )
}

const$env <- ros

CJS_tcovAge.out <- nimbleMCMC(code = CJS_tcovAge,
                           data = data, constants = const,inits = initsFunction,
                           monitors = c('phi.mu', 'p.mu','ran.phi', 'ran.p','slope', 'phi.sd', 'p.sd'),
                           summary=T,samplesAsCodaMCMC = T,WAIC = F,
                           nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)



ggout <-  CJS_tcovAge.out$samples %>% ggs() %>%
  filter(!grepl('ran.',Parameter))

ggout %>% filter(Parameter==c('phi.mu[1]','phi.mu[2]') )%>%
  ggs_traceplot()

ggout %>% filter(Parameter=='phi.sd') %>%
  ggs_traceplot()


# M-array -----------------------------------------------------------------

# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
    nind <- dim(CH)[1]
    n.occasions <- dim(CH)[2]
    m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
    # Calculate the number of released individuals at each time period
    for (t in 1:n.occasions){
        m.array[t,1] <- sum(CH[,t])
    }
    for (i in 1:nind){
        pos <- which(CH[i,]!=0)
        g <- length(pos)
        for (z in 1:(g-1)){
            m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
        } #z
    } #i
    # Calculate the number of individuals that is never recaptured
    for (t in 1:n.occasions){
        m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
    }
    out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
    return(out)
}



# 7.10.2. Time-dependent models
# Specify model in BUGS language
CJS_MArray <- nimbleCode({
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
        phi[t] ~ dunif(0, 1)         # Priors for survival
        p[t] ~ dunif(0, 1)           # Priors for recapture
    }
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
        marr[t,1:n.occasions] ~ dmulti(pr[t,1:n.occasions], rel[t])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
        q[t] <- 1-p[t]                # Probability of non-recapture
        pr[t,t] <- phi[t]*p[t]
        # Above main diagonal
        for (j in (t+1):(n.occasions-1)){
            pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
        } #j
        # Below main diagonal
        for (j in 1:(t-1)){
            pr[t,j] <- 0
        } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
        pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } #t
})

# Create the m-array from the capture-histories
marr <- marray(obs)

# Bundle data
dataM <- list(marr = marr)
constM <- list(n.occasions = dim(marr)[2], rel = rowSums(marr))

# Initial values
initsFunction <- function(){list(phi = runif(dim(marr)[2]-1, 0, 1), p = runif(dim(marr)[2]-1, 0, 1))}

start <- Sys.time()
CJS_MArray.out <- nimbleMCMC(code = CJS_MArray,
                              data = dataM, constants = constM,inits = initsFunction,
                              monitors = c("phi", "p"  ),
                              summary=T,samplesAsCodaMCMC = T,WAIC = F,
                             nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)
CJS_MArray.dur <- Sys.time()-start

CJS_MArray.out$samples %>% ggs %>% ggs_caterpillar(family = 'phi',sort = F)


tmp2 <- CJS_tt.out$samples %>% ggs() %>%
    filter(Parameter %in% paste0('phi[',1:25,']'))
tmp3 <- CJS_MArray.out$samples %>% ggs() %>%
    filter(Parameter %in% paste0('phi[',1:25,']'))
ggs_caterpillar(list(tt=tmp2,Marray=tmp3),sort=F,horizontal = F)+
    coord_cartesian(x=c(0.25,1))
