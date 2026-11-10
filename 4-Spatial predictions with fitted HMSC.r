# Spatial predictions with fitted HMSC models (presence/absence or abundance)

rm(list = ls())
source("./0-Load libraries.R")

# This script generates spatial predictions from a fitted HMSC model for either
# presence/absence ("mpa") or abundance ("abu"). It saves both full posterior
# distributions and posterior means per year.

# ------------------------------------------------------------------------------
# General settings
# ------------------------------------------------------------------------------

model_type <- "mpa"  # "mpa" for presence/absence, "abu" for abundance
thin <- 4000
nSamples <- 250
nChains <- 4
SpatRes <- "0.05"
Depthrange <- "0_1000"

# Output directories
dir_save_posteriors <- if (model_type == "mpa") {
  "./Outputs/Presence_Absence_posteriors/"
} else {
  "./Outputs/Abundance_posteriors/"
}

dir_save_mean_posterior <- if (model_type == "mpa") {
  "./Outputs/Presence_Absence_mean_posterior/"
} else {
  "./Outputs/Abundance_mean_posterior/"
}

dir.create(dir_save_posteriors, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_save_mean_posterior, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load fitted model
# ------------------------------------------------------------------------------

model_file <- paste0(
  "./Models/model_", model_type, "_thin_", thin,
  "_samples_", nSamples, "_chains_", nChains, ".RData"
)
load(model_file)

# ------------------------------------------------------------------------------
# Load environmental grid data
# ------------------------------------------------------------------------------

grid_env <- get(load(
  paste0("./Inputs_HMSC/Hindcast_Env_grid_rect", SpatRes, "_WMED_", Depthrange, ".RData")
))
# Drop potential duplicate rows per (X, Y, YEAR) to avoid redundant computation
grid_env <- dplyr::distinct(grid_env, X, Y, YEAR, .keep_all = TRUE)
grid_env$Depth <- (-1) * grid_env$Depth

# ------------------------------------------------------------------------------
# Spatial predictions by year
# ------------------------------------------------------------------------------

# Note: This loop can be parallelized by year on HPC. It can be slow and storage-heavy.
# We recommend running it on a High-Performance Computing (HPC) cluster with sufficient resources. Generated ouptuts (>500 GB) for both presence/absence and abundance models can be requested from the authors.
vec_year <- 1999:2021

for (i in seq_along(vec_year)) {
  cat("i:", i, "\n")
  grid_env_year <- grid_env[grid_env$YEAR == vec_year[i], ]

  t1 <- Sys.time()

  xy.grid <- as.matrix(cbind(grid_env_year$X, grid_env_year$Y))
  XData.grid <- data.frame(
    Depth = grid_env_year$Depth,
    SST = grid_env_year$SST,
    SBT = grid_env_year$botTemp,
    SSS = grid_env_year$so,
    Chla = grid_env_year$chl,
    Fishing = grid_env_year$FPI_tot,
    Gravity = grid_env_year$Gravity,
    stringsAsFactors = TRUE
  )

  Gradient <- prepareGradient(m, XDataNew = XData.grid, sDataNew = list(grid.cell = xy.grid))

  # Posterior predictive distribution
  predY <- predict(
    m, Gradient = Gradient, predictEtaMean = TRUE, expected = TRUE,
    mcmcStep = 1, nParallel = 4, useSocket = TRUE
  )

  predY_array <- array(
    unlist(predY),
    dim = c(nrow(predY[[1]]), ncol(predY[[1]]), length(predY))
  )
  saveRDS(predY_array, file = paste0(dir_save_posteriors, model_type, "_posteriors_", vec_year[i], ".rds"))

  rm(predY_array)
  gc()

  # Posterior mean prediction
  EpredY <- Reduce("+", predY) / length(predY)
  saveRDS(EpredY, file = paste0(dir_save_mean_posterior, model_type, "_mean_posterior_", vec_year[i], ".rds"))

  rm(list = c("predY", "EpredY"))
  gc()

  t2 <- Sys.time()
  print(t2 - t1)
}