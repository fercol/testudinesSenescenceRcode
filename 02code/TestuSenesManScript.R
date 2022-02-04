# ============================= CODE METADATA ================================ #
# AUTHOR: Fernando Colchero
# SCRIPT TYPE: Main code
# DATE CREATED: 2022-02-01
# DESCRIPTION: Demo code for reproducibility of manuscript Da Silva et al.
# ============================== START CODE ================================== #
# ======================== #
# ==== GENERAL SETUP: ====
# ======================== #
# Setup working directory (point it to the parent folder 
# "testudinesSenescenceRcode":
# setwd("Directory where the parent folder is")

# List of installed libraries:
instPacks <- installed.packages()[, 1]

# Install missing libraries:
if (!"snowfall" %in% instPacks) {
  install.packages("snowfall")
}

if (!"RColorBrewer" %in% instPacks) {
  install.packages("RColorBrewer")
}

if (!"phytools" %in% instPacks) {
  install.packages("phytools")
}

if (!"ape" %in% instPacks) {
  install.packages("ape")
}

if (!"caper" %in% instPacks) {
  install.packages("caper")
}

if (!"BaSTA.ZIMS" %in% instPacks) {
  install.packages("05packages/BaSTA.ZIMS_1.0.2.tar.gz", 
                   type = "source", repos = NULL)
}

if (!"paramDemo" %in% instPacks) {
  install.packages("05packages/paramDemo_1.0.0.tar.gz", 
                   type = "source", repos = NULL)
}

# Load libraries:
library(snowfall)
library(phytools)
library(ape)
library(caper)
library(nlme)
library(BaSTA.ZIMS)
library(paramDemo)

# Source functions:
source("02code/TestuSenesfunctions.R")

# ==================== #
# ==== LOAD DATA: ====
# ==================== #
# Load phylogentic tree:
phyloAll <- read.newick("03data/phylo/1-s2.0-S1055790316304316-mmc2.nwk.txt")

# Load covariate data:
covarsDat <- read.table("03data/tables/covariateData.txt", 
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Load response variables:
responseDat <- read.table("03data/tables/responseVars.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)

# Example of data for BaSTA:
bastaDat <- read.table("03data/tables/testBaSTAData.txt", sep = "\t", 
                       header = TRUE)

# ========================= #
# ==== BaSTA ANALYSIS: ====
# ========================= #
# Run BaSTA on test data for Siler model:
out <- bastaZIMS(bastaDat, model = "GO", shape = "bathtub", ncpus = 4, 
                 nsim = 4, niter = 21001, burnin = 1000, parallel = TRUE)

# Check goodness of fit:
plot(out, plot.type = "gof")

# Extract Siler mortality parameters:
theta <- out$coefficients[, 1]

# Extract summary statistics for life expectancy:
leSumar <- c(Mean = mean(out$PS$nocov$Ex), SD = sd(out$PS$nocov$Ex),
             Lower = quantile(out$PS$nocov$Ex, 0.025, names = FALSE),
             Upper = quantile(out$PS$nocov$Ex, 0.975, names = FALSE))

# Run parallel routine to extract aging rate summary statistics:
arSumar <- CalcARsumStats(out, ncpus = 4, sxlevs = c(0.5, 0.2), 
                          model = "GO", shape = "bathtub")

# =================================== #
# ==== PHYLOGENETIC REGRESSIONS: ====
# =================================== #
# ----------------------------------------------------- #
# Example of regression for life expectancy in females:
# ----------------------------------------------------- #
# Extract point estimate for life expectancy and aging rates:
resDat <- responseDat[, c("species", "Mean.EX.Female", "Mean.AR.Female")]

# Rename columns:
colnames(resDat) <- c("species", "ex", "ar")

# Variable transformation for life expectancy:
resDat$ex <- log(resDat$ex)

# Prepare regression data:
regrDat <- PrepRegrData(responseData = resDat, 
                         predictorData = covarsDat, 
                         predictors = c('habitat', 'hibernation', 
                                        'reprFemales', 'afr', 
                                        "weightFemales", "bodyWeightDiff",
                                        "reprOutDiff"),
                         phyloAll = phyloAll)


# Run regression:
outEX1F <- RunBayesPGLS(formula = ex ~ habitat + hibernation +
                          reprFemales + afr + log(weightFemales) + 
                          bodyWeightDiff + reprOutDiff, 
                        data = regrDat$data, Sigma = regrDat$Sigma, 
                        nsim = 4, ncpus = 4)
# Verify output:
summary(outEX1F)

# Plot posterior densities for parameters:
plot(outEX1F, plot.type = 'density')

# Plot diagnostics:
plot(outEX1F, plot.type = 'diagnostics')


# ----------------------------------------------------- #
# Example of regression for life expectancy in males:
# ----------------------------------------------------- #
# Extract point estimate for life expectancy and aging rates:
resDat <- responseDat[, c("species", "Mean.EX.Male", "Mean.AR.Male")]

# Rename columns:
colnames(resDat) <- c("species", "ex", "ar")

# Variable transformation for life expectancy:
resDat$ex <- log(resDat$ex)

# Prepare regression data:
regrDat <- PrepRegrData(responseData = resDat, 
                        predictorData = covarsDat, 
                        predictors = c('habitat', 'hibernation', 
                                       'reprMales', 'afr', 
                                       "weightMales", "bodyWeightDiff",
                                       "reprOutDiff"),
                        phyloAll = phyloAll)


# Run regression:
outEX1M <- RunBayesPGLS(formula = ex ~ habitat + hibernation +
                          reprMales + afr + log(weightMales) + 
                          bodyWeightDiff + reprOutDiff, 
                        data = regrDat$data, Sigma = regrDat$Sigma, 
                        nsim = 4, ncpus = 4)
# Verify output:
summary(outEX1M)

# Plot posterior densities for parameters:
plot(outEX1M, plot.type = 'density')

# Plot diagnostics:
plot(outEX1M, plot.type = 'diagnostics')

