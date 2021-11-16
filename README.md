# Estimer la survie et la densité grâce à un suivi d'individus marqués: une introduction avec une approche Bayesienne

**durée:** 4 h

**nb max de participant:** 25

## préalables: 
- utiliser R
- bonne base en stats (glm binomial, modèle mixte)
- connaissance minimale en bayésien


## Plan du workshop

1. Introduction: rapide mise à niveau sur le bayésien et l'utilisation de nimble
2. Introduction: c'est quoi du Capture-Marquage-Recapture 
3. Estimer la survie: le modèle Cormak-Jolly-Seber
4. Estimer la survie et +: les modèles multi-state
5. Estimer la densité: CMR en population fermée 
6. Vers la fin et plus loin encore 


## pour le workshop

Avoir un laptop

Installer R et RStudio

Installer nimble:
https://r-nimble.org/html_manual/cha-installing-nimble.html

Installer les packages suivants:
`install.packages(c("tidyverse", "ggmcmc", "here", "nimbleEcology"))`

S'assurer que le code suivant roule sans erreur:
```
      library(nimble)
      code <- nimbleCode({
      y ~ dnorm(0,1)
      })
      model <- nimbleModel(code)
      cModel <- compileNimble(model)
```
