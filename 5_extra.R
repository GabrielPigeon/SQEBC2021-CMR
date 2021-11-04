# using nimble ecology for CJS
library(nimbleEcology)

CJS_tt_eco <- nimbleCode({
    for (t in 1:(n.occasions-1)){
        p[t] ~ dunif(0, 1)         # Priors for time-spec. recapture
        phi[t] ~ dunif(0, 1)        # Priors for time-spec. survival
    }
    p[n.occasions] ~ dunif(0, 1)
    # Likelihood
    for (i in 1:nind){
        obs[i,f[i]:n.occasions] ~ dCJS_vv(phi[f[i]:(n.occasions-1)], p[f[i]:n.occasions], tl[i])
    } #i
})

initsFunction <- function(){
    list(phi = runif(ncol(obs)-1,0.7,0.9),
         p = runif(ncol(obs),.4,.7)
    )
}

dataEco <- list(obs=cbind(obs,0))
constEco <- list(n.occasions=ncol(dataEco$obs),
                 nind=nrow(obs),
                 f=yrbirth,
                 tl = ncol(dataEco$obs)-yrbirth+1)

start <- Sys.time()
CJS_tt_eco.out <- nimbleMCMC(code = CJS_tt_eco,
                             data = dataEco, constants = constEco,
                             monitors = c("phi", "p"),
                             summary=T,samplesAsCodaMCMC = T,WAIC = F,
                             nburnin = myburn,thin = mythin, niter = mynit, nchains = 3)
CJS_tt_eco.dur <- Sys.time()-start


tmp1 <- CJS_tt_eco.out$samples %>% ggs() %>%
    filter(Parameter %in% paste0('phi[',1:25,']'))


# multi-site example ------------------------------------------------------

multisite <- nimbleCode({

    # -------------------------------------------------
    # Parameters:
    # phiA: survival probability site A
    # phiB: survival probability site B
    # phiC: survival probability site B
    # psiAA: movement probability from site A to site A (reference)
    # psiAB = psiA[1]: movement probability from site A to site B
    # psiAC = psiA[2]: movement probability from site A to site C
    # psiBA = psiB[1]: movement probability from site B to site A
    # psiBB: movement probability from site B to site B (reference)
    # psiBC = psiB[2]: movement probability from site B to site C
    # psiCA = psiC[1]: movement probability from site C to site A
    # psiCB = psiC[2]: movement probability from site C to site B
    # psiCC: movement probability from site C to site C (reference)
    # pA: recapture probability site A
    # pB: recapture probability site B
    # pC: recapture probability site C
    # -------------------------------------------------
    # States (z):
    # 1 alive at A
    # 2 alive at B
    # 2 alive at C
    # 3 dead
    # Observations (y):
    # 1 not seen
    # 2 seen at A
    # 3 seen at B
    # 3 seen at C
    # -------------------------------------------------

    # Priors
    phiA ~ dunif(0, 1)
    phiB ~ dunif(0, 1)
    phiC ~ dunif(0, 1)
    pA ~ dunif(0, 1)
    pB ~ dunif(0, 1)
    pC ~ dunif(0, 1)
    # transitions: multinomial logit
    # normal priors on logit of all but one transition probs
    for (i in 1:2){
        lpsiA[i] ~ dnorm(0, sd = 1000)
        lpsiB[i] ~ dnorm(0, sd = 1000)
        lpsiC[i] ~ dnorm(0, sd = 1000)
    }
    # constrain the transitions such that their sum is < 1
    for (i in 1:2){
        psiA[i] <- exp(lpsiA[i]) / (1 + exp(lpsiA[1]) + exp(lpsiA[2]))
        psiB[i] <- exp(lpsiB[i]) / (1 + exp(lpsiB[1]) + exp(lpsiB[2]))
        psiC[i] <- exp(lpsiC[i]) / (1 + exp(lpsiC[1]) + exp(lpsiC[2]))
    }
    # last transition probability  (because it must sum to 1)
    psiA[3] <- 1 - psiA[1] - psiA[2]
    psiB[3] <- 1 - psiB[1] - psiB[2]
    psiC[3] <- 1 - psiC[1] - psiC[2]

    # alternative using Dirichlet priors (add alpha = c(1, 1, 1) to constants given to model)
    # psiA[1:3] ~ ddirch(alpha[1:3])
    # psiB[1:3] ~ ddirch(alpha[1:3])
    # psiC[1:3] ~ ddirch(alpha[1:3])

    # probabilities of state z(t+1) given z(t)
    gamma[1,1] <- phiA * psiA[1]
    gamma[1,2] <- phiA * psiA[2]
    gamma[1,3] <- phiA * psiA[3]
    gamma[1,4] <- 1 - phiA
    gamma[2,1] <- phiB * psiB[1]
    gamma[2,2] <- phiB * psiB[2]
    gamma[2,3] <- phiB * psiB[3]
    gamma[2,4] <- 1 - phiB
    gamma[3,1] <- phiC * psiC[1]
    gamma[3,2] <- phiC * psiC[2]
    gamma[3,3] <- phiC * psiC[3]
    gamma[3,4] <- 1 - phiC
    gamma[4,1] <- 0
    gamma[4,2] <- 0
    gamma[4,3] <- 0
    gamma[4,4] <- 1

    # probabilities of y(t) given z(t)
    omega[1,1] <- 1 - pA     # Pr(alive A t -> non-detected t)
    omega[1,2] <- pA         # Pr(alive A t -> detected A t)
    omega[1,3] <- 0          # Pr(alive A t -> detected B t)
    omega[1,4] <- 0          # Pr(alive A t -> detected C t)
    omega[2,1] <- 1 - pB     # Pr(alive B t -> non-detected t)
    omega[2,2] <- 0          # Pr(alive B t -> detected A t)
    omega[2,3] <- pB         # Pr(alive B t -> detected B t)
    omega[2,4] <- 0          # Pr(alive B t -> detected C t)
    omega[3,1] <- 1 - pC     # Pr(alive C t -> non-detected t)
    omega[3,2] <- 0          # Pr(alive C t -> detected A t)
    omega[3,3] <- 0          # Pr(alive C t -> detected B t)
    omega[3,4] <- pC         # Pr(alive C t -> detected C t)
    omega[4,1] <- 1          # Pr(dead t -> non-detected t)
    omega[4,2] <- 0          # Pr(dead t -> detected A t)
    omega[4,3] <- 0          # Pr(dead t -> detected B t)
    omega[4,4] <- 0          # Pr(dead t -> detected C t)

    # likelihood
    for (i in 1:N){
        # latent state at first capture
        z[i,first[i]] <- y[i,first[i]] - 1
        for (t in (first[i]+1):K){
            # z(t) given z(t-1)
            z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
            # y(t) given z(t)
            y[i,t] ~ dcat(omega[z[i,t],1:4])
        }
    }
})

