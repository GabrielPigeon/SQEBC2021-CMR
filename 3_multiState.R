library(tidyverse)
library(ggmcmc)
library(nimble)
library(coda)
library(here)
load(here('data','reindeer.Rdata'))
mythin=1
myburn=500
mynit=myburn+mythin*1000



##   mvt between sites ------------------

MS_model1 <- nimbleCode({
    # Parameters:
    # phiA: survival probability at site A
    # phiB: survival probability at site B
    # psiAB: movement probability from site A to site B
    # psiBA: movement probability from site B to site A
    # pA: recapture probability at site A
    # pB: recapture probability at site B
    # States (S):
    # 1 alive at A
    # 2 alive at B
    # 3 dead
    # Observations (O):
    # 1 seen at A
    # 2 seen at B
    # 3 not seen

    # Priors and constraints
    for (t in 1:(n.occasions-1)){
        phiA[t] <- mean.phi[1]
        phiB[t] <- mean.phi[2]
        psiAB[t] <- mean.psi[1]
        psiBA[t] <- mean.psi[2]
        pA[t] <- mean.p[1]
        pB[t] <- mean.p[2]
    }
    for (u in 1:2){
        mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
        mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
        mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
    }

    # Define state-transition and observation matrices
    for (i in 1:nind){
        # Define probabilities of state S(t+1) given S(t)
        for (t in f[i]:(n.occasions-1)){
            ps[1,i,t,1] <- phiA[t] * (1-psiAB[t])
            ps[1,i,t,2] <- phiA[t] * psiAB[t]
            ps[1,i,t,3] <- 1-phiA[t]
            ps[2,i,t,1] <- phiB[t] * psiBA[t]
            ps[2,i,t,2] <- phiB[t] * (1-psiBA[t])
            ps[2,i,t,3] <- 1-phiB[t]
            ps[3,i,t,1] <- 0
            ps[3,i,t,2] <- 0
            ps[3,i,t,3] <- 1

            # Define probabilities of O(t) given S(t)
            po[1,i,t,1] <- pA[t]
            po[1,i,t,2] <- 0
            po[1,i,t,3] <- 1-pA[t]
            po[2,i,t,1] <- 0
            po[2,i,t,2] <- pB[t]
            po[2,i,t,3] <- 1-pB[t]
            po[3,i,t,1] <- 0
            po[3,i,t,2] <- 0
            po[3,i,t,3] <- 1
        } #t
    } #i

    # Likelihood
    for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- y[i,f[i]]
        for (t in (f[i]+1):n.occasions){
            # State process: draw S(t) given S(t-1)
            z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,1:3])
            # Observation process: draw O(t) given S(t)
            y[i,t] ~ dcat(po[z[i,t], i, t-1,1:3])
        } #t
    } #i
})

# Function to create known latent states z
known.state.ms <- function(ms, notseen){
# notseen: label for ‘not seen’
    state <- ms
    state[state==notseen] <- NA
    for (i in 1:dim(ms)[1]){
m <- min(which(!is.na(state[i,])))
state[i,m] <- NA
}
return(state)
}

# Function to create initial values for unknown z
ms.init.z <- function(ch, f){
for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
states <- max(ch, na.rm = TRUE)
known.states <- 1:(states-1)
v <- which(ch==states)
ch[-v] <- NA
ch[v] <- sample(known.states, length(v), replace = TRUE)
return(ch)
}

# Bundle data
data <- list(y = MS1_obs,
                  z = known.state.ms(MS1_obs, 3))
const=list(f = yrbirth,
           n.occasions = dim(MS1_obs)[2], nind = dim(MS1_obs)[1])

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1),
                         mean.psi = c(max(rnorm(1, .3, .2),0),rnorm(1, .5, .2)),
                         mean.psi = runif(2, 0, 1),
                         mean.p = runif(2, 0, 1),
                         z = ms.init.z(MS1_obs, yrbirth))}

start <- Sys.time()
ms1.out <- nimbleMCMC(code = MS_model1,
                         data = data, constants = const,inits = inits,
                         monitors = c("mean.phi", "mean.psi", "mean.p"),
                         summary=T,samplesAsCodaMCMC = T,
                      nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)
ms1.dur <- Sys.time()-start

ms1.out$summary$all.chains
ms1.out$samples %>% ggs() %>%ggs_traceplot(family = 'psi')


# repro state -------------------------------------------------------------

MS_model2 <- nimbleCode({
    # States (S):
    # 1 alive  Juv
    # 2 alive breeder
    # 3 alive-nonbreeder
    # Observations (O):
    # 1 seen nocalf
    # 2 seen calf
    # 3 not seen

    # Priors and constraints
    for (t in 1:(n.occasions-1)){
        ranef.p[t]~dnorm(0,sd=sd.p)
        logit(pj[t]) <- mu.p[1]+ranef.p[t]
        logit(pb[t]) <- mu.p[2]+ranef.p[t]
        logit(pnb[t]) <- mu.p[3]+ranef.p[t]
        ranef.nb[t]~dnorm(0,sd.nb)
        logit(psiBNb[t]) <- mu.psiBNb+ranef.nb[t]
    }

    sd.p~dunif(0,5)
    sd.nb~dunif(0,5)
    mu.psiBNb~ dunif(0, 1)
    psiNbB~ dunif(0, 1)

    for (u in 1:3){
        mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
        mu.p[u] ~ dlogis(0, 1)      # Priors for mean state-spec. recapture
        mean.p[u] <- ilogit(mu.p[u])
    }

    psiJ[1] <- 0
    for(a in 2:nbAge){
        psiJ[a] ~ dunif(0, 1)    # Priors for transitions
    }


    # Define state-transition and observation matrices
    for (i in 1:nind){
        # Define probabilities of state S(t+1) given S(t)
        for (t in f[i]:(n.occasions-1)){
            ps[1,i,t,1] <- mean.phi[1] * (1-psiJ[age[i,t]])
            ps[1,i,t,2] <- mean.phi[1] * psiJ[age[i,t]]
            ps[1,i,t,3] <- 0
            ps[1,i,t,4] <- 1-mean.phi[1]

            ps[2,i,t,1] <- 0
            ps[2,i,t,2] <- mean.phi[2] * (1-psiBNb[t])
            ps[2,i,t,3] <- mean.phi[2] *psiBNb[t]
            ps[2,i,t,4] <- 1-mean.phi[2]

            ps[3,i,t,1] <- 0
            ps[3,i,t,2] <- mean.phi[3]*psiNbB
            ps[3,i,t,3] <- mean.phi[3]*(1-psiNbB)
            ps[3,i,t,4] <- 1-mean.phi[3]


            ps[4,i,t,1] <- 0
            ps[4,i,t,2] <- 0
            ps[4,i,t,3] <- 0
            ps[4,i,t,4] <- 1

            # Define probabilities of O(t) given S(t)
            po[1,i,t,1] <- pj[t]
            po[1,i,t,2] <- 0
            po[1,i,t,3] <- 1-pj[t]

            po[2,i,t,1] <- 0
            po[2,i,t,2] <- pb[t]
            po[2,i,t,3] <- 1-pb[t]

            po[3,i,t,1] <- pnb[t]
            po[3,i,t,2] <- 0
            po[3,i,t,3] <- 1-pnb[t]

            po[4,i,t,1] <- 0
            po[4,i,t,2] <- 0
            po[4,i,t,3] <- 1
        } #t
    } #i

    # Likelihood
    for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
            # State process: draw S(t) given S(t-1)
            z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,1:4])
            # Observation process: draw O(t) given S(t)
            y[i,t] ~ dcat(po[z[i,t], i, t-1,1:3])
        } #t
    } #i

    # derided
    afr[1] <- psiJ[1]
    for(a in 2:nbAge){
    afr[a] <- prod(1-psiJ[1:(a-1)])*psiJ[a]
    }
    afr[8] <- prod(1-psiJ[1:(7)])
})


# Bundle data
ageC <-c(1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7,7,
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7)
know.state <- ifelse(MS2_obs==2,2,NA)
for(i in 1:nrow(know.state)){
    know.state[i,yrbirth[i]] <- 1
    if(MS2_obs[i,yrbirth[i]+1] == 1) know.state[i,yrbirth[i]+1] <- 1
    }

data <- list(y = MS2_obs,
             z = know.state)

const=list(f = yrbirth,
           age=apply(age, 1:2, function(x) ageC[x]),
           n.occasions = dim(MS2_obs)[2],
           nind = dim(MS2_obs)[1]
)
const$nbAge=max(const$age,na.rm = T)

# Initial values
inits <- function(){
    l <- list(mean.phi = runif(3, 0, 1),
              psiJ = runif(max(ageC), 0, 1),
              psiNbB=runif(1,0,1),mu.psiBNb=runif(1,0,1),
            mu.p = rnorm(3, 0, 1))
    l$z <- ifelse(is.na(know.state),2,NA)
return(l)
}

start <- Sys.time()
ms2.out <- nimbleMCMC(code = MS_model2,
                      data = data, constants = const,inits = inits,
                      monitors = c("mean.phi", "psiBNb",'psiNbB','psiBNb','mu.psiBNb','psiJ', "mean.p",'afr'),
                      summary=T,samplesAsCodaMCMC = T,
                      nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)
ms2.dur <- Sys.time()-start

ms2.out$summary$all.chains
ms2.out$samples %>% ggs() %>%ggs_traceplot(family = 'phi')
ms2.out$samples %>% ggs() %>%ggs_traceplot(family = 'p\\[')
ms2.out$samples %>% ggs() %>%ggs_traceplot(family = 'psiBNb\\[3\\]')
ms2.out$samples %>% ggs() %>%ggs_traceplot(family = 'psiNbB')

ms2.out$samples %>% ggs() %>%ggs_caterpillar(family = 'psiJ')

ms2.out$samples %>% ggs() %>%
    ggs_caterpillar(family = 'afr',sort = F,horizontal = F)+
    scale_y_discrete(name='Age de la premiere repro',labels=c(1:8))+
    theme(axis.text.x = element_text(angle=0))+
    labs(x='Probabilité')



# carcasse recovery -------------------------------------------------------


MS3_rec <- nimbleCode({
    # Priors and constraints
    mean.p ~ dlogis(0, 1)     # Prior for mean survival
    mean.r ~ dlogis(0, 1)     # Prior for mean recovery
    sd.p ~dunif(0,10)
    sd.r ~dunif(0,10)
    Omega[1:nbAge, 1:nbAge]  ~ dwish(R[1:nbAge, 1:nbAge], nbAge+1)

    for(a in 1:nbAge){
        mu.s[a] ~ dlogis(0, 1)     # Prior for mean recapture
        mean.s[a] <- ilogit(mu.s[a])
        xi[a] ~dunif(0,10)
    }

    for (t in 1:(n.occasions-1)){
        ranef.r[t]~dnorm(0,sd=sd.r)
        logit(r[t]) <- mean.r+ranef.r[t]
        ranef.p[t]~dnorm(0,sd=sd.p)
        logit(p[t]) <- mean.p+ranef.p[t]

        eps.yr[t,1:nbAge]~ dmnorm(zero[1:nbAge], Omega[1:nbAge, 1:nbAge])
        for(a in 1:nbAge){
            ranef.yr[t,a] <-  xi[a]*eps.yr[t,a]
        }
    }

    for(i in 1:nind){
        for (t in f[i]:(n.occasions-1)){
            logit(s[i,t]) <- mu.s[age[i,t]]+ranef.yr[t,age[i,t]]
        }
    }

    # Define state-transition and observation matrices
    for (i in 1:nind){
        # Define probabilities of state S(t+1) given S(t)
        for (t in f[i]:(n.occasions-1)){
            ps[1,i,t,1] <- s[i,t]#*F[t]
            ps[1,i,t,2] <- (1-s[i,t])*r[t]
            ps[1,i,t,3] <- (1-s[i,t])*(1-r[t])
            ps[2,i,t,1] <- 0
            ps[2,i,t,2] <- 0
            ps[2,i,t,3] <- 1
            ps[3,i,t,1] <- 0
            ps[3,i,t,2] <- 0
            ps[3,i,t,3] <- 1

            # Define probabilities of O(t) given S(t)
            po[1,i,t,1] <- p[t]
            po[1,i,t,2] <- 0
            po[1,i,t,3] <- 1-p[t]
            po[2,i,t,1] <- 0
            po[2,i,t,2] <- 1
            po[2,i,t,3] <- 0
            po[3,i,t,1] <- 0
            po[3,i,t,2] <- 0
            po[3,i,t,3] <- 1
        } #t
    } #i

    # Likelihood
    for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- y[i,f[i]]
        for (t in (f[i]+1):n.occasions){
            # State process: draw S(t) given S(t-1)
            z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,1:3])
            # Observation process: draw O(t) given S(t)
            y[i,t] ~ dcat(po[z[i,t], i, t-1,1:3])
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

# Bundle data
ageC <-c(1, 2, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5,
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)
data <- list(y = MS3_obs)

const=list(f = yrbirth,
           age=apply(age, 1:2, function(x) ageC[x]),
           n.occasions = dim(MS3_obs)[2],
           nind = dim(MS3_obs)[1]
)

const$nbAge=max(const$age,na.rm = T)
const$zero=rep(0,const$nbAge)
const$R = diag(const$nbAge)

# Initial values

ld.init2 <- function(ch, f){
    v2 <- which(ch==2, arr.ind = T)
    ch[v2] <- 2
    for (i in 1:nrow(v2)){
        ifelse(v2[i,2]!=ncol(ch), ch[v2[i,1], (v2[i,2]+1):ncol(ch)] <- 3, next)}
    for (i in 1:nrow(ch)){
        m <- max(which(ch[i,]==1))
        ch[i,f[i]:m] <- 1
    }
    for (i in 1:nrow(v2)){
        u1 <- min(which(ch[v2[i,1],]==1))
        ch[v2[i],u1:(v2[i,2]-1)] <- 1
    }
    for (i in 1:nrow(ch)){
        for (j in f[i]:ncol(ch)){
            if(is.na(ch[i,j])==1) ch[i,j] <- 1
        }
        ch[i,f[i]] <- NA
    }
    return(ch)
}



inits <- function(){list(mean.p = runif(1, 0, 1), mean.f = runif(1, 0, 1),
                         mu.s = rnorm(6, 0, 1), mean.r = runif(1, 0, 1),
                         sd.p=runif(1,0,1),sd.r=runif(1,0,1),
                         ranef.p=rnorm(ncol(MS3_obs)-1,0,0.5),
                         ranef.r=rnorm(ncol(MS3_obs)-1,0,0.5),
                         xi=runif(6,0,1),
                         eps.yr=matrix(rnorm((ncol(MS3_obs)-1)*6,0,1),nrow =ncol(MS3_obs)-1,ncol=6 ),
                         Omega=inverse(diag(runif(6,0,1))),
                         z = ld.init2(MS3_obs, yrbirth))}



start <- Sys.time()
ms3.out <- nimbleMCMC(code = MS3_rec,
                      data = data, constants = const,inits = inits(),
                      monitors = c('mean.p','mean.r','mu.s', #'Omega', 'Sigma.raw',
                                   'rho','Sigma','ranef.yr','sd.p','sd.r'),
                      summary=T,samplesAsCodaMCMC = T,
                      nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)
ms3.dur <- Sys.time()-start

ms3.out$summary$all.chains
ms3.out$samples %>% ggs() %>%ggs_traceplot(family = 'mu')
ms3.out$samples %>% ggs() %>%ggs_traceplot(family = 'mean')
ms3.out$samples %>% ggs() %>%ggs_traceplot(family = 'Sigma')
