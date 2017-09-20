#############################
### Construct data input ----
#############################
setwd("~/Dropbox/WorkPersonal_transfers/PhD/Kautikeino/Compendium")
load("./Data/richness_all.rda")
library(dplyr)

# isolate relevant env variables
richness_all_JDSM <- richness_all %>%
  filter(!is.na(area),
         !is.na(elevation),
         !is.na(shoreline_complex))

# edit env variables
Xpre <- richness_all_JDSM %>%
  transmute(log_area = log(area),elev = elevation,shore = shoreline_complex)
Xpre$intercept <- c(rep(1,327))
Xpre <- Xpre %>%
  dplyr::select(4,1,2,3)

# Run loop on glms to determine coefficients for fixed effects of env variables
coefs1 <- data.frame(Intercept=numeric(0),Area=numeric(0),Elevation=numeric(0),Shoreline=numeric(0))
species_names <- c("Burbot","ArcticCharr",
  "BrownTrout","Grayling","Whitefish","Perch","Pike")
for (i in species_names)
{
  coefs1 <- rbind(coefs1,c("1","1","1","1"))
}
colnames(coefs1) <- c("intercept","area","elevation","shore")
rownames(coefs1) <- species_names

initial_val <- data.frame(Intercept=numeric(0),Area=numeric(0),Elevation=numeric(0),Shoreline=numeric(0))
for (i in species_names)
{
  JSDM_example <- glm(data=richness_all,family=binomial,richness_all[,i] ~ 
                        scale(log(area)) + scale(eurolst_bio01) + shoreline_complex)
  initial_val <- rbind(initial_val,coefficients(JSDM_example))
}
colnames(initial_val) <- c("intercept","area","elevation","shore")
rownames(initial_val) <- species_names




#################
### Run JSDM ----
#################

packages.needed <-
  setdiff(
    c('R2jags', 'MASS', 'MCMCpack', 'abind', 'random', 'mclust'),
    rownames(installed.packages())
  )
if(length(packages.needed)) install.packages(packages.needed)

X <- data.matrix(Xpre)[,-1]
coefs <- data.matrix(coefs1)[,-1]
Occur <- model.matrix(~.+0,dplyr::select(richness_all_JDSM,c(2:8)))[,2:8]

n.chains <- 5
n.iter <- 150000
n.burn <- 11000
n.thin <- 40

df <- 2

model_name <- 'new_model'
source('./R/fit_JSDM.R')

save(new_model,file="new_model.rda")
load("new_model.rda")

#############################
#### Save JDSM results ------
#############################

SumEnv1
SumEnv <- SUMMARY(Beta,mean)
rownames(SumEnv) <- colnames(Occur)
colnames(SumEnv) <- c("Intercept",colnames(coefs))
SumEnvSD <- SUMMARY(Beta,sd)
rownames(SumEnvSD) <- colnames(Occur)
colnames(SumEnvSD) <- c("Intercept",colnames(coefs))

SumFish1
SumFish <- SUMMARY(EnvRho,mean)
colnames(SumFish) <- colnames(Occur)
rownames(SumFish) <- colnames(Occur)
SumFishSD <- SUMMARY(EnvRho,sd)
colnames(SumFishSD) <- colnames(Occur)
rownames(SumFishSD) <- colnames(Occur)


save(SumEnv,file="./Data/SumEnv.rda")
save(SumEnvSD,file="./Data/SumEnvSD.rda")
save(SumFish,file="./Data/SumFish.rda")
save(SumFishSD,file="./Data/SumFishSD.rda")

load("./Data/SumEnv.rda")
load("./Data/SumFish.rda")


Diagnose(Beta, 'effn')

