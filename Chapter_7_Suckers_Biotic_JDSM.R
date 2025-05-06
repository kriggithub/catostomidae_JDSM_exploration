#### K. Riggin May 6, 2025
#### Joint Species Distribution Modeling Biotic Interactions



####Step 1: Load in packages and data

library(tidyverse)
library(Hmsc)
library(corrplot)


data <- read.csv("suckerspresencedata.csv")


# separate model construction and validation sites
data_C <- data %>% filter(val_cond == "C")
data_V <- data %>% filter(val_cond == "V")



# make sure to select either all common names or all scientific names



Y <- as.matrix(data_C %>% select(golden.redhorse:torrent.sucker))


# selecting aggregated eco region as variable
XData <- data.frame(aggecoreg = as.factor(data_C$AG_ECO9), temp = data_C$ws_temp, erosion = data_C$erod_ws)

studyDesign <- data.frame(site_id = as.factor(data_C$SITE_ID))
rL <- HmscRandomLevel(units = studyDesign$site_id)

XFormula <- ~ aggecoreg + temp + erosion
XFormula2 <- ~ 1



model <- Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit", studyDesign = studyDesign, ranLevels = list(site_id = rL), YScale = T)
model2 <- Hmsc(Y = Y, XData = XData, XFormula = XFormula2, distr = "probit", studyDesign = studyDesign, ranLevels = list(site_id = rL), YScale = T)



nChains = 2
samples = 2000
thin = 10
transient = 5000
verbose = 1000
nParallel = 2
modelsample <- sampleMcmc(model, thin = thin, samples = samples, transient = transient,
                          nChains = nChains, verbose = verbose, nParallel = nParallel, 
                          initPar = "fixed effects")


modelsample2 <- sampleMcmc(model2, thin = thin, samples = samples, transient = transient,
                          nChains = nChains, verbose = verbose, nParallel = nParallel, 
                          initPar = "fixed effects")


modelpost <- convertToCodaObject(modelsample, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
# for beta parameters
plot(modelpost$Beta) # good convergence for some, terrible for others
effectiveSize(modelpost$Beta) # ESS high for some, very low for ohers
gelman.diag(modelpost$Beta, multivariate = FALSE)$psrf # relatively close to 1 for some, far for others
summary(modelpost$Beta)





preds <- computePredictedValues(modelsample, expected = F)
modelfit <- evaluateModelFit(hM = modelsample, predY = preds)
modelfit




VP <- computeVariancePartitioning(modelsample)
VP


postBeta <- getPostEstimate(modelsample, parName="Beta")

plotBeta(modelsample, post=postBeta, param="Support", supportLevel = 0.95, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))



OmegaCor <- computeAssociations(modelsample)
supportLevel <- 0.95
toPlot <- ((OmegaCor[[1]]$support > supportLevel) +
             (OmegaCor[[1]]$support < (1-supportLevel)) > 0) * OmegaCor[[1]]$mean
corrplot(toPlot, method = "color", col = c("grey", "white", "black"))



OmegaCor2 <- computeAssociations(modelsample2)
supportLevel2 <- 0.95
toPlot2 <- ((OmegaCor2[[1]]$support > supportLevel2) +
             (OmegaCor2[[1]]$support < (1-supportLevel2)) > 0) * OmegaCor2[[1]]$mean
corrplot(toPlot2, method = "color", col = c("grey", "white", "black"))



biPlot(modelsample, etaPost = getPostEstimate(modelsample, "Eta"),
       lambdaPost = getPostEstimate(modelsample, "Lambda"), colVar = 2)


biPlot(modelsample2, etaPost = getPostEstimate(modelsample2, "Eta"),
       lambdaPost = getPostEstimate(modelsample2, "Lambda"), colVar = 2)





