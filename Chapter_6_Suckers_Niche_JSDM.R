#### K. Riggin May 5, 2025
#### Joint Species Distribution Modeling Traits and Phylogeny


# https://www.youtube.com/watch?v=u07eFE3Uqtg&t=1445s&ab_channel=EarthLabCUBoulder
# https://github.com/admahood/ff_study 

####Step 1: Load in packages and data
library(ape)
library(tidyverse)
library(Hmsc)


data <- read.csv("suckerspresencedata.csv")
tree <- read.nexus("mito_tree.nex")
traits <- read.csv("trait_subset.csv")
com_to_sci <- read.csv("com_&_sci_names.csv")





# select data that includes species from traits and phylogeny
matchedcomnames <- as.character(traits$COMMONNAME)
matchedscinames <- as.character(traits$SCINAME)
data <- data %>% select(SITE_ID:val_cond, any_of(matchedcomnames))



name_vector <- setNames(com_to_sci$Species, com_to_sci$COMMONNAME)

data <- data %>%
  rename_with(~ name_vector[.x], .cols = any_of(names(name_vector)))



# separate model construction and validation sites
data_C <- data %>% filter(val_cond == "C")
data_V <- data %>% filter(val_cond == "V")



# make sure to select either all common names or all scientific names



Y <- as.matrix(data_C %>% select(carpiodes.carpio:moxostoma.rupiscartes))
XData <- data.frame(temp = data_C$ws_temp)
TrData <- data.frame(temptol = traits$MAXTEMP)
rownames(TrData) <- matchedscinames #rownames of TrData must match species names

studyDesign <- data.frame(site_id = as.factor(data_C$SITE_ID))
rL <- HmscRandomLevel(units = studyDesign$site_id)

XFormula <- ~ temp
TrFormula <- ~ temptol


model <- Hmsc(Y = Y, XData = XData, XFormula = XFormula, TrData = TrData, TrFormula = TrFormula,
              phyloTree = tree, distr = "probit", studyDesign = studyDesign, ranLevels = list(site_id = rL))



nChains = 2
samples = 1000
thin = 50
transient = 5000
verbose = 1000
nParallel = 2
modelsample <- sampleMcmc(model, thin = thin, samples = samples, transient = transient,
                          nChains = nChains, verbose = verbose, nParallel = nParallel,
                          initPar = "fixed effects")


modelpost <- convertToCodaObject(modelsample, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
# for beta parameters
plot(modelpost$Beta) # very poor convergence, only shows every other speciess
effectiveSize(modelpost$Beta) # ESS very low
gelman.diag(modelpost$Beta, multivariate = FALSE)$psrf # bad values, far from one except for hypentelium.nigricans
summary(modelpost$Beta)


preds <- computePredictedValues(modelsample, expected = F)
evaluateModelFit(hM = modelsample, predY = preds)



postGamma = getPostEstimate(modelsample, parName="Gamma")
plotGamma(modelsample, post=postGamma, param="Support", supportLevel = 0.25, covNamesNumbers = c(T,F), trNamesNumbers = c(T,F), colorLevels = 3)



Gradient <- constructGradient(modelsample, focalVariable = "temp")
predY <- predict(modelsample, Gradient = Gradient, expected = T)
plotGradient(modelsample, Gradient = Gradient, pred = predY, measure = "S",
             showData = T)
plotGradient(modelsample, Gradient = Gradient, pred = predY, measure = "T", index = 2, showData = T)



VP <- computeVariancePartitioning(modelsample)
VP





####Validate model predictions
# setup predictions for validation dataset
studyDesign2 <- data.frame(site_id = as.factor(data_V$SITE_ID))
rL2 <- HmscRandomLevel(units = studyDesign$site_id)
XData.grid <- data.frame(temp = data_V$ws_temp)
Gradient.grid <- prepareGradient(modelsample, XDataNew = XData.grid,
                                 sDataNew = list(site_id = rL))
predYnr <- predict(modelsample, Gradient = Gradient.grid)



# posterior occurrence probability for each site
EpredY <- as.data.frame(apply(abind(predYnr, along = 3), c(1,2), mean))
data_V$observed <- rowSums(data_V[, 17:29])
EpredY$SITE_ID <- data_V$SITE_ID
posteriordata <- data_V %>% select(SITE_ID, ws_temp, observed) %>% left_join(EpredY, by = "SITE_ID")
posteriordata$expected <- rowSums(posteriordata[, 3:15])
posteriordata$oe <- posteriordata$observed/posteriordata$expected

mean(posteriordata$oe)


ggplot(posteriordata, aes(y = oe)) +
  geom_boxplot()

hist(posteriordata$oe)



# give a sense of what we were missing from each trait from and also phylogeny (what trait and what percent of species were missing)
# give the taxonomy a try instead of phylogeny


# of the 28 species in original data set, only 23 had data on temperature tolerance or substrate/current preference, and of those, 3 did not have temperature tolerance data


# phylogeny only had 13 of the original 28 species from the dataset, looking back seems like there are small speeling differences like ii vs i for white suckers (Caotsomtus commersonii)


# limiting factor might be phylogeny


# given the issue of subsetting everything
# predict joint probability (probability for each species)
# sum predicted occurence probabilities across species which will end up being distribution of E at that site
# could just take the median value of E and then divide it by O




