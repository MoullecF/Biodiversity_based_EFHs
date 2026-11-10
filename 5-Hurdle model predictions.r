# Hurdle model predictions (presence/absence × abundance)

rm(list = ls())
gc()

# Load one fitted model to get species names
load("./Models/model_abu_thin_4000_samples_250_chains_4.RData")
sp.names <- colnames(m$Y)

# Number of posterior samples
nSamples <- 1000

# Confidence interval width (approx. 95% by default)
ci <- function(x, z = 1.96) {
  mean(x, na.rm = TRUE) + z * sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))) -
    (mean(x, na.rm = TRUE) - z * sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
}

# ------------------------------------------------------------------------------
# Inputs and templates
# ------------------------------------------------------------------------------

list.PA <- list.files(path = "./Outputs/Presence_Absence_posteriors", full.names = TRUE)
list.AB <- list.files(path = "./Outputs/Abundance_posteriors", full.names = TRUE)

template <- readRDS(list.PA[1])
preds_t <- array(NA, dim = c(dim(template), length(list.PA)))

# ------------------------------------------------------------------------------
# Hurdle predictions by year
# ------------------------------------------------------------------------------

vec_year <- 1999:2021

for (i in seq_along(vec_year)) {

  year <- vec_year[i]
  cat("Year:", year, "\n")
  t1 <- Sys.time()

  Post_PA <- readRDS(list.PA[grep(pattern = as.character(year), list.PA)])
  Post_AB <- readRDS(list.AB[grep(pattern = as.character(year), list.AB)])

  for (post in 1:nSamples) {
    cat("Post:", post, fill = TRUE)
    # Multiply PA by abundance (undo log with exp)
    preds_t[, , post] <- Post_PA[, , post] * exp(Post_AB[, , post])
  }

  # Point estimates (mean, sd, CI width)
  preds_t_mean <- apply(preds_t, c(1, 2), mean, na.rm = TRUE)
  preds_t_sd <- apply(preds_t, c(1, 2), sd, na.rm = TRUE)
  preds_t_ci <- apply(preds_t, c(1:2), ci)

  saveRDS(preds_t_mean, paste0("./Outputs/Hurdle_prediction/Mean/Mean_Pred_AB_Cond_Pres_", year, ".rds"))
  saveRDS(preds_t_sd, paste0("./Outputs/Hurdle_prediction/SD/SD_Pred_AB_Cond_Pres_", year, ".rds"))
  saveRDS(preds_t_ci, paste0("./Outputs/Hurdle_prediction/Ci/Ci_Pred_AB_Cond_Pres_", year, ".rds"))

  t2 <- Sys.time()
  print(t2 - t1)
  gc()
}

# Generated outputs represent more than 1GB of data. They can be requested from the authors.

# ------------------------------------------------------------------------------
# Spatial grid and XY coordinates
# ------------------------------------------------------------------------------

grid_env <- get(load("./Inputs_HMSC/Hindcast_Env_grid_rect0.05_WMED_0_1000.RData"))
# Drop potential duplicate rows per (X, Y, YEAR) to avoid redundant computation
grid_env <- dplyr::distinct(grid_env, X, Y, YEAR, .keep_all = TRUE)
xy <- grid_env[grid_env$YEAR == 1999, c("X", "Y")]

WD_SpatialPred <- "./Outputs/Hurdle_prediction/Mean/"

# ------------------------------------------------------------------------------
# Build spatio-temporal raster stack from yearly predictions
# ------------------------------------------------------------------------------

list.Pred <- list.files(WD_SpatialPred, full.names = TRUE)

# Initialize stacks
r.stack.year <- stack()

for (i in 1:length(list.Pred)) {

  Pred_year_prov <- readRDS(list.Pred[i])
  Pred_year_prov <- cbind(Pred_year_prov, xy)
  Pred_year_prov <- relocate(Pred_year_prov, X, Y)
  year <- str_sub(list.Pred[i], start = 57, end = 60)

  r.stack <- stack()

  cat("Year:", year, "\n")

  for (j in 3:dim(Pred_year_prov)[2]) {

    XYZ_data_sp <- Pred_year_prov[, c(1:2, j)]
    sp_name <- colnames(XYZ_data_sp)[3]
    r <- rasterFromXYZ(
      xyz = XYZ_data_sp,
      crs = "+proj=longlat +datum=WGS84",
      res = c(0.05, 0.05)
    )
    names(r) <- paste0(sp_name, "_", year)
    r.stack <- stack(r.stack, r)
  }

  r.stack.year <- stack(r.stack.year, r.stack)
}

# save(r.stack.year, file = "./Outputs/Spatio_temporal_prediction/r_stack_Hurdle_0.05_0_1000_19992021.Rdata")
