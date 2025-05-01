#### K. Riggin April 22, 2025
#### Single Species Distribution Modeling
#### Chosen species: White Sucker (present in 199 sites with max count of 169)


####Step 1: Load in packages, data, and preliminary inspection
library(Hmsc)
library(mapview)
library(abind)
library(tidyverse)
library(ggspatial)
library(glmmTMB)

data <- read.csv("suckerspresencedata.csv")
data$elevation <- data$elevation/100

# separate model construction and validation sites
data_C <- data %>% filter(val_cond == "C")
data_V <- data %>% filter(val_cond == "V")


# Look at location of sites
mapview(data, xcol = "long", ycol = "lat", crs = 4269, grid = FALSE)
mapview(data_C, xcol = "long", ycol = "lat", crs = 4269, grid = FALSE)
mapview(data_V, xcol = "long", ycol = "lat", crs = 4269, grid = FALSE)



# preliminary relationships between abundance, occurrence and covariates
  # elevation
ggplot(data_C, aes(x = elevation, y = WHITE.SUCKER)) +
  geom_point() +
  geom_smooth(method = "lm")
elevationoccurance <- glm(WHITE.SUCKER ~ elevation, data = data_C, family = binomial(link = 'logit'))
summary(elevationoccurance)

  # precipitation
ggplot(data_C, aes(x = ws_precip, y = WHITE.SUCKER)) +
  geom_point() +
  geom_smooth(method = "lm")
precipoccurance <- glm(WHITE.SUCKER ~ ws_precip, data = data_C, family = binomial(link = 'logit'))
summary(precipoccurance)

  # temperature
ggplot(data_C, aes(x = ws_temp, y = WHITE.SUCKER)) +
  geom_point() +
  geom_smooth(method = "lm")
tempoccurance <- glm(WHITE.SUCKER ~ ws_temp, data = data_C, family = binomial(link = 'logit'))
summary(tempoccurance)

  # runoff
ggplot(data_C, aes(x = ws_runoff, y = WHITE.SUCKER)) +
  geom_point() +
  geom_smooth(method = "lm")
runoffoccurance <- glm(WHITE.SUCKER ~ ws_runoff, data = data_C, family = binomial(link = 'logit'))
summary(runoffoccurance)

  # slope
ggplot(data_C, aes(x = slope, y = WHITE.SUCKER)) +
  geom_point() +
  geom_smooth(method = "lm")
slopeoccurance <- glm(WHITE.SUCKER ~ slope, data = data_C, family = binomial(link = 'logit'))
summary(slopeoccurance)

  # erosion
ggplot(data_C, aes(x = erod_ws, y = WHITE.SUCKER)) +
  geom_point() +
  geom_smooth(method = "lm")
erodoccurance <- glm(WHITE.SUCKER ~ erod_ws, data = data_C, family = binomial(link = 'logit'))
summary(erodoccurance)




####Step 2: Subset Data for Modeling
# select temperature and erosion to have relevance to trait data (max temp/substrate preference)

XData <- data.frame(temp = data_C$ws_temp, erosion = data_C$erod_ws)

# select White Suckers to model (Catostomus commersonii)
Y <- as.matrix(data_C$WHITE.SUCKER)
colnames(Y) <- "Catostomus_commersonii"
# bind together coordinates to include spatial random effect
latlong <- as.matrix(cbind(data_C$lat, data_C$long))


# spatial graph for initial inspection
ggplot(data, aes(x = long, y = lat, fill = ws_temp)) +
  borders("state", colour = "gray80", fill = "gray90") +
  geom_point(shape = 21, color = "black", size = 2) +
  scale_fill_gradient(low = "blue", high = "red") +
  coord_fixed(xlim = c(-95, -65), ylim = c(32, 50))


ggplot(data, aes(x = long, y = lat, fill = erod_ws)) +
  borders("state", colour = "gray80", fill = "gray90") +
  geom_point(shape = 21, color = "black", size = 2) +
  scale_fill_gradient(low = "blue", high = "red") +
  coord_fixed(xlim = c(-95, -65), ylim = c(32, 50))



####Step 3: Setup model and run MCMC chain
studyDesign <- data.frame(site_id = as.factor(data_C$SITE_ID))
rownames(latlong) <- studyDesign[,1]
rL <- HmscRandomLevel(sData = latlong)
XFormula <- ~ temp + erosion

# use probit model for occurance data
model <- Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
              studyDesign = studyDesign, ranLevels = list(site_id = rL))


# MCMC parameters
nChains = 2
samples = 1000
thin = 10
transient = 5000
verbose = 1000
# MCMC sampling with first fitting models in maximum-likelihood framework
modelsample <- sampleMcmc(model, thin = thin, samples = samples, transient = transient,
                           nChains = nChains, verbose = verbose, initPar = "fixed effects")



####Step 4: Evaluate MCMC convergence
modelpost <- convertToCodaObject(modelsample, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F)) # refer to species and covars by names instead of numbers
  # for beta parameters
plot(modelpost$Beta) # convergence looks good
effectiveSize(modelpost$Beta) # ESS relatively high
gelman.diag(modelpost$Beta, multivariate = FALSE)$psrf # relatively close to 1
summary(modelpost$Beta)
  # for spatial scale parameter
plot(modelpost$Alpha[[1]]) # convergence is very poor
effectiveSize(modelpost$Alpha[[1]]) 
gelman.diag(modelpost$Alpha[[1]], multivariate = FALSE)$psrf 







####Step 5: Evaluate explanatory & predictive power of model
  # Explanatory power
preds <- computePredictedValues(modelsample, expected = F)
evaluateModelFit(hM = modelsample, predY = preds)
# AUC = 0.915
# TjurR2 = 0.515
# Indicates relatively high explanatory power for determining presences and absences




# Construct gradient plots for Beta Parameters for model predictions
par(mfrow  = c(1,2))
tempGradient <- constructGradient(modelsample, focalVariable = "temp")
temppred <- predict(modelsample, Gradient = tempGradient, expected = T)
plotGradient(modelsample, tempGradient, pred = temppred, measure = "Y", index = 1, showData = T)

eroGradient <- constructGradient(modelsample, focalVariable = "erosion")
eropred <- predict(modelsample, Gradient = eroGradient, expected = T)
plotGradient(modelsample, eroGradient, pred = eropred, measure = "Y", index = 1, showData = T)






# # cross validation ########## Should have not included validation references? Explanatory measurements went down significantly
# partition <- createPartition(model.sample, nfolds = 2, column = "site_id")
# predscv <- computePredictedValues(model.sample, partition = partition)
# evaluateModelFit(hM = model.sample, predY = predscv)





# Predict distribution in validation sites
latlong.grid <- as.matrix(cbind(data_V$lat, data_V$long))
XData.grid <- data.frame(temp = data_V$ws_temp, erosion = data_V$erod_ws)
Gradient.grid <- prepareGradient(modelsample, XDataNew = XData.grid,
                                 sDataNew = list(site_id = latlong.grid))
predYnr <- predict(modelsample, Gradient = Gradient.grid)


# posterior occurrence probability for each site
EpredY <- as.data.frame(apply(abind(predYnr, along = 3), c(1,2), mean))
EpredY$SITE_ID <- data_V$SITE_ID
posteriordata <- data_V %>% select(SITE_ID, lat, long, ws_temp, erod_ws, WHITE.SUCKER) %>% left_join(EpredY, by = "SITE_ID")
posteriordata$WHITE.SUCKER <- ifelse(posteriordata$WHITE.SUCKER == 1, "present", "absent")


# summary(posteriordata$long)
# summary(posteriordata$lat)
ggplot(posteriordata, aes(x = long, y = lat, fill = `Catostomus commersonii`)) +
  borders("state", colour = "gray80", fill = "gray90") +
  geom_point(shape = 21, color = "black", size = 2) +
  scale_fill_gradient(low = "blue", high = "red") +
  coord_fixed(xlim = c(-90, -67), ylim = c(33, 50))

ggplot(posteriordata, aes(x = long, y = lat, fill = WHITE.SUCKER)) +
  borders("state", colour = "gray80", fill = "gray90") +
  geom_point(shape = 21, color = "black", size = 2) +
  coord_fixed(xlim = c(-90, -67), ylim = c(33, 50))














