#### K. Riggin April 22, 2025
#### Single Species Distribution Modeling
#### Chosen species: White Sucker (present in 199 sites with max count of 169)


####Step 1: Load in packages and datasets from cleaning
library(Hmsc)
library(mapview)
library(abind)
env_dat <- read.csv("environmental_table_subset.csv")
fish_dat <- read.csv("fish_count_subset.csv")
reference_dat <- read.csv("reference_table_subset.csv")

data <- env_dat %>%
  left_join(reference_dat %>% select(SITE_ID, AG_ECO9, ref_cond, val_cond), by = "SITE_ID") %>% 
  left_join(fish_dat, by = "SITE_ID")
  
# seperate SAP and NAP to visualize
data_SAP <- data %>%
  filter(AG_ECO9 == "SAP")
data_NAP <- data %>%
  filter(AG_ECO9 == "NAP")


# separate reference sites
data_ref <- data %>% 
  filter(ref_cond == "REF")

data_ref_whitesuckers <- data_ref %>% 
  filter(WHITE.SUCKER != 0)


data_nonref <- data %>% 
  filter(ref_cond == "NON-REF")


# Look at location of sites
mapview(data, xcol = "long", ycol = "lat", crs = 4269, grid = FALSE)
mapview(data_SAP, xcol = "long", ycol = "lat", crs = 4269, grid = FALSE)
mapview(data_NAP, xcol = "long", ycol = "lat", crs = 4269, grid = FALSE)
mapview(data_ref, xcol = "long", ycol = "lat", crs = 4269, grid = FALSE)
mapview(data_ref_whitesuckers, xcol = "long", ycol = "lat", crs = 4269, grid = FALSE)




####Step 2: Subset Data for Modeling
# select elevation and precipitation as environmental covariates
XData <- data.frame(elevation = data_ref$elevation, precip = data_ref$ws_precip) 
# select White Suckers to model (Catostomus commersonii)
Y <- as.matrix(data_ref$WHITE.SUCKER)
colnames(Y) <- "Catostomus commersonii"
# bind together coordinates
latlong <- as.matrix(cbind(data_ref$lat, data_ref$long))



####Step 3: Setup model
studyDesign <- data.frame(site_id = as.factor(data_ref$SITE_ID))
rownames(latlong) <- studyDesign[,1]
rL <- HmscRandomLevel(sData = latlong)
XFormula <- ~ elevation + precip

# lognormal poisson for count data
model <- Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
              studyDesign = studyDesign, ranLevels = list(site_id = rL))




nChains = 2
samples = 1000
thin = 10
transient = 5000
verbose = 1000

model.sample <- sampleMcmc(model, thin = thin, samples = samples, transient = transient,
                           nChains = nChains, verbose = verbose, initPar = "fixed effects")



####Step 4: Evaluate MCMC convergence (increased thinning and transient size due to bad initial convergence)
model.post <- convertToCodaObject(model.sample)
summary(model.post$Beta)
effectiveSize(model.post$Beta)
gelman.diag(model.post$Beta, multivariate = FALSE)$psrf
plot(model.post$Beta)



####Step 5: evaluate model fit
preds <- computePredictedValues(model.sample, expected = F)
evaluateModelFit(hM = model.sample, predY = preds)



# variance partitioning
groupnames <- c("elevation", "precipitation")
group <- c(1,2)
VP <- list()
VP <- computeVariancePartitioning(model.sample, group = group, groupnames = groupnames)
VP



# see how parameters translate into model predictions
par(mfrow  = c(1,2))
Gradient <- constructGradient(model.sample, focalVariable = "precip")
predY <- predict(model.sample, Gradient = Gradient, expected = T)
plotGradient(model.sample, Gradient, pred = predY, measure = "Y", index = 1, showData = T)


Gradient2 <- constructGradient(model.sample, focalVariable = "elevation")
predY2 <- predict(model.sample, Gradient = Gradient2, expected = T)
plotGradient(model.sample, Gradient2, pred = predY2, measure = "Y", index = 1, showData = T)






# cross validation ########## Should have not included validation references? Explanatory measurements went down significantly
partition <- createPartition(model.sample, nfolds = 2, column = "site_id")
predscv <- computePredictedValues(model.sample, partition = partition)
evaluateModelFit(hM = model.sample, predY = predscv)





# predict distribution in non-reference sites
latlong.grid <- as.matrix(cbind(data_nonref$lat, data_nonref$long))
XData.grid <- data.frame(elevation = data_nonref$elevation, precip = data_nonref$ws_precip)
Gradient.grid <- prepareGradient(model.sample, XDataNew = XData.grid,
                                 sDataNew = list(site_id = latlong.grid))
predYnr <- predict(model.sample, Gradient = Gradient.grid)
# predictive mean value over posterior
EpredY <- apply(abind(predYnr, along = 3), c(1,2), mean)
# posterior probability that the count is non-zero
EpredO <- apply(abind(predYnr, along = 3), c(1,2), 
                FUN = function(a) {mean(a > 0)})
plotGradient(model.sample, Gradient.grid, predY = predYnr, measure = "Y", index = 1, showData = T)




