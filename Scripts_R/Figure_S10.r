###############################################################################
# Title: Figure S10 - Gamma Diagnostics (ESS + PSRF)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

source("Scripts_R/0-Load libraries.R")

# -----------------------------------------------------------------------------
# Load models and convert to coda objects
# -----------------------------------------------------------------------------

mPA <- get(load("Models/model_mpa_thin_4000_samples_250_chains_4.RData"))
mABU <- get(load("Models/model_abu_thin_4000_samples_250_chains_4.RData"))

mpost_PA <- convertToCodaObject(mPA)
mpost_ABU <- convertToCodaObject(mABU)

# -----------------------------------------------------------------------------
# Plot: ESS and PSRF for gamma
# -----------------------------------------------------------------------------

png(filename = file.path("Figures", "Figure_S10.png"), width = 10, height = 10, units = "in", res = 400)
par(mfrow = c(2, 2))

hist(
  effectiveSize(mpost_PA$Gamma),
  main = "ess(gamma) - PA model",
  xlab = "Effective sample size"
)

vioplot(
  gelman.diag(mpost_PA$Gamma, multivariate = FALSE)$psrf,
  main = "psrf(gamma) - PA model"
)

hist(
  effectiveSize(mpost_ABU$Gamma),
  main = "ess(gamma) - ABU model",
  xlab = "Effective sample size"
)

vioplot(
  gelman.diag(mpost_ABU$Gamma, multivariate = FALSE)$psrf,
  main = "psrf(gamma) - ABU model"
)

dev.off()
