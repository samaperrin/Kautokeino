##################################
### Import and configure data ----
##################################

setwd("~/Dropbox/WorkPersonal_transfers/PhD/Kautikeino/Compendium")
load("./Data/richness_all.rda")
load("./Data/distance_matrix.rda")
source("./R/db_connect.R")
library(HMSC)


# Start Auto matrix - distance between sites
KautoDist <- data.frame(cmdscale(distance_matrix),rownames(distance_matrix));colnames(distance_matrix)<- c("X1","X2","LocationID")

## Create X matrix - environmental variables

# Filter down variables to those we have data for on everything
richness_all_HMSC <- richness_all %>%
  filter(!is.na(area),
         !is.na(elevation),
         !is.na(shoreline_complex),
         locationID %in% rownames(distance_matrix))
richness_all_HMSC <- richness_all_HMSC[!duplicated(richness_all_HMSC$locationID), ]

# Isolate relevant variables and get into correct format
X <- richness_all_HMSC %>%
  mutate(log_area = log(area),
         elevation = log(elevation),temp=log(eurolst_bio10)) %>%
         dplyr::select(log_area,elevation,shoreline_complex,temp)
rownames(X) <- seq(1,317,1)
dimnames(X) <- list(seq(1,317,1), c("log_area","elevation","shoreline_complex","temperature"))
# Make extra one for prediction matrices later on
X_WithTemp <- X

### Finish Auto matrix to match locations coordinates up
KautoDist <- KautoDist[match(richness_all_HMSC$locationID, KautoDist$rownames.distance_matrix.),]
rownames(KautoDist) <- seq(1,317,1)
KautoDist <- KautoDist %>%
  dplyr::select(3,1,2) %>%
  rename(sampling_unit = rownames.distance_matrix.)

# Create Y matrix - community species composition
# Isolate occurence variables and make them integers
Y <- richness_all_HMSC %>%
  dplyr::select(c(2:8))
for(i in 1:ncol(Y)){
  Y[,i] <- as.integer(as.character(Y[,i]))
}
YMat <- as.matrix(Y)

# Create Tr matrix - traits
# Need to do this manually, based off info from Per-Arne/Anders
Tr <- t(as.data.frame(cbind(c(1,1,1,0,1,1,1),
                            c(0,1,1,1,0,1,0),
                            c(0,1,0,1,1,1,0))))
rownames(Tr) <- c("Littoral","Profundal","Pelagic")
dimnames(Tr) <- list(c("Littoral","Profundal","Pelagic"),seq(1,7,1))




# Put in HMSC format
KautoHMSCData <- as.HMSCdata(Y = YMat, X = X_WithTemp, Tr = Tr,scaleTr=FALSE,Auto = KautoDist)

# Create priors and define parameters
KautoPrior <- as.HMSCprior(KautoHMSCData)
KautoParam <- as.HMSCparam(KautoHMSCData, KautoPrior)

# Run model (takes at least 6 hours)
model <- hmsc(KautoHMSCData_WithTemp2, family = "probit", niter = 20000, nburn = 10000,
              thin = 10)

# or just load it from below
# load("./Data/kauto_hmsc2.rda")

##############################################
### Data analysis - Environmental effects ----
##############################################

# Analysis of posterior distribution

# Trace and density plots
# NB: Need large margin figures
mixing <- as.mcmc(KautoModel_WithTemp2, parameters = "paramX")
par(mfrow = c(3,2))
par(mar = c(2, 2, 2, 2))
plot(mixing)

### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing)

### Draw beanplots
library(beanplot)
par(mfrow = c(1,1))
par(mar = c(13, 4, 2, 2))
beanplot(mixingDF, las = 2)

### True values
truth <- as.vector(KautoParam_WithTemp2$param$paramX)

### Average
average <- apply(KautoModel_WithTemp2$results$estimation$paramX, 1:2, mean)

### 95% confidence intervals
CI.025 <- apply(KautoModel_WithTemp2$results$estimation$paramX, 1:2, quantile,
                probs = 0.025)
CI.975 <- apply(KautoModel_WithTemp2$results$estimation$paramX, 1:2, quantile,
                probs = 0.975)
CI <- cbind(as.vector(CI.025), as.vector(CI.975))

### Draw confidence interval plots
par(mar = c(4, 4, 2, 2))
plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI, truth), type = "n",
     xlab = "", ylab = "", main="paramX")
abline(h = 0,col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2],
       code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15, cex = 1.5)
points(1:nrow(CI), truth, col = "red", pch = 19)

### Summary table
paramXCITable <- cbind(unlist(as.data.frame(CI.025)),
                       unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("2.5%", "Mean", "97.5%")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")

##################################################################################
### Additional output from a model with traits and phylogenetic correlations #####
##################################################################################

### Mixing object for paramTr
mixing <- as.mcmc(KautoModel_WithTemp2, parameters = "paramTr")

### Draw trace and density plots for all combination of parameters
par(mfrow = c(3,2))
plot(mixing)

### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing)

### Draw boxplot for each parameters
par(mfrow = c(1,1))
par(mar = c(12, 4, 2, 2))
boxplot(mixingDF, las = 2)

### Draw beanplots
library(beanplot)
par(mar = c(10, 4, 2, 2))
beanplot(mixingDF, las = 2)

### True values
truth <- as.vector(KautoParam_WithTemp2$param$paramTr)

### Average
average <- apply(KautoModel_WithTemp2$results$estimation$paramTr, 1:2, mean)

### 95% confidence intervals
CI.025 <- apply(KautoModel_WithTemp2$results$estimation$paramTr, 1:2, quantile,
                probs = 0.025)
CI.975 <- apply(KautoModel_WithTemp2$results$estimation$paramTr, 1:2, quantile,
                probs = 0.975)
CI <- cbind(as.vector(CI.025), as.vector(CI.975))

### Draw confidence interval plots
plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI, truth), type = "n",
     xlab = "", ylab = "", main = "paramX")
abline(h = 0, col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2],
       code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15, cex = 1.5)
points(1:nrow(CI), truth, col = "red", pch = 19)

### Summary table
paramXCITable <- cbind(unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.025)),
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")




###############################################
### Data analysis - Variation partitioning ----
###############################################

variationPartKauto_WithTemp2 <- variPart(KautoModel_WithTemp2, 
              c(rep("Area",2),"Elevation","Shoreline","Temperature"))
save(variationPartKauto_WithTemp2,file="variationPartKauto_WithTemp2.rda")
par(mar=c(6,3,5,1))
barplot(t(variationPartKauto_WithTemp2$variPart), las=2, cex.names=0.75, cex.axis=0.75,
        legend.text=paste(colnames(variationPartKauto_WithTemp2$variPart)," ",
        signif(100*colMeans(variationPartKauto_WithTemp2$variPart),2),"%",sep=""),
        args.legend=list(y=1.2, xjust=1, horiz=F, bty="n",cex=0.75))
###########################################
### Data analysis - Correlation matrix ----
###########################################

### Draw estimated correlation matrix
library(corrplot)
corMat <- corRandomEff(KautoModel_WithTemp2, cor = TRUE)
averageCor <- apply(corMat[, , , 1], 1:2, mean)
for (i in 1:7) {
  averageCor[i,i] <- 0
}
maxAvg <- max(c(max(averageCor),-min(averageCor)))
averageCor <- averageCor*floor(1/(maxAvg))
corrplot(averageCor, method = "color",
         col = colorRampPalette(c("red", "white", "blue"))(200))

### Draw chord diagram
library(circlize)
corMat <- corRandomEff(KautoModel_WithTemp2, cor = TRUE)
colMat <- matrix(NA, nrow = nrow(averageCor), ncol = ncol(averageCor))
colMat[which(averageCor > 0.6, arr.ind = TRUE)] <- "blue"
colMat[which(averageCor < -0.6, arr.ind = TRUE)] <- "red"
chordDiagram(averageCor, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey",col=colMat)



Ymean_withtemp2 <- apply(KautoModel_WithTemp2$data$Y, 2, mean)
R2_withtemp2 <- Rsquared(KautoModel_WithTemp2, averageSp=FALSE)
plot(Ymean_withtemp2, R2_withtemp2, pch=19, main=paste('Mean R2 value over species',
                                   signif(mean(R2_withtemp2),2)))


########## Predicting community composition with 2.5 degree rise in temperature

# Create version of matrix where temperature is not scaled
X_NoScale <- richness_all_HMSC %>%
  mutate(log_area = log(area),
         elevation = log(elevation),temp=eurolst_bio10)
X_NoScale$intercept <- rep(1,317)
X_NoScale <- X_NoScale %>%
  dplyr::select(intercept,log_area,elevation,shoreline_complex,temp)


### Prepare simulation of system where temperature higher
meanX <- apply(X_NoScale, 2, mean)
meanX[1] <- 0
sc <- apply(KautoHMSCData$X, 2, sd) / apply(X_NoScale, 2, sd)
sc[1] <- 1

### Create environmental matrix with increased temperature
X_IncTemp <- richness_all_HMSC %>%
  mutate(log_area = scale(log(area)),
         elevation = scale(log(elevation)),temp=scale(log(eurolst_bio10))+1) %>%
  dplyr::select(log_area,elevation,shoreline_complex,temp)

# New trait matrix for this purpose
X_IncTemp$Intercept <- rep(1,317)
X_IncTemp_Mat <- as.matrix(dplyr::select(X_IncTemp,5,1,2,3,4))
Intercept <- rep(1,7)
Tr_pred <- rbind(Intercept,Tr)

n <- 900

### Predictions

predVecR<-array(NA,dim=c(nrow(X_IncTemp_Mat),ncol(KautoModel_WithTemp2$data$Y),n))

set.seed(3)
samp <- sample(1:dim(KautoModel_WithTemp2$results$estimation$paramX)[3],
               n,replace=TRUE)
for (i in 1:n) {
  i2 <- samp[i]
  
  ## Define objects
  paramXModel <- KautoModel_WithTemp2$results$estimation$paramX[,,i2]
  paramTrModel <- KautoModel_WithTemp2$results$estimation$paramTr[,,i2]
  latentAutoModel <- matrix(KautoModel_WithTemp2$results$estimation$latentAuto[i2,], 
                        ncol = ncol(KautoModel_WithTemp2$results$estimation$latentAuto))
  paramLatentAutoModel <- matrix(KautoModel_WithTemp2$results$estimation$paramLatentAuto[i2,], 
                             ncol = ncol(KautoModel_WithTemp2$results$estimation$paramLatentAuto))
  
  ### Simulate community
  predVecR[,,i] <- communitySimul(X=as.matrix(X_IncTemp_Mat),paramX = paramXModel,
                    Tr=Tr_pred,paramTr = paramTrModel,Auto=KautoDist)$probMat
}

as.data.frame(rbind(colnames(Y),round(apply(predVecR,2,mean))))
apply(YMat,2,mean)



##########################################################################################
##########################################################################################
##########################################################################################


