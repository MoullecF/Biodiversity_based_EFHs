###############################################################################
# Title: Figure S4 - Beta Diagnostics (ESS + PSRF)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

source("Scripts_R/0-Load libraries.R")

# -----------------------------------------------------------------------------
# Load models and convert to coda objects
# -----------------------------------------------------------------------------

mPA <- get(load("../Models/model_mpa_thin_4000_samples_250_chains_4.RData"))
mABU <- get(load("../Models/model_abu_thin_4000_samples_250_chains_4.RData"))

mpost_PA <- convertToCodaObject(mPA)
mpost_ABU <- convertToCodaObject(mABU)

# -----------------------------------------------------------------------------
# Plot: ESS and PSRF for beta
# -----------------------------------------------------------------------------

png(filename = file.path("Figures", "Figure_S4.png"), width = 10, height = 10, units = "in", res = 400)
par(mfrow = c(2, 2))

hist(
  effectiveSize(mpost_PA$Beta),
  main = "ess(beta) - PA model",
  xlab = "Effective sample size"
)

vioplot(
  gelman.diag(mpost_PA$Beta, multivariate = FALSE)$psrf,
  main = "psrf(beta) - PA model"
)

hist(
  effectiveSize(mpost_ABU$Beta),
  main = "ess(beta) - ABU model",
  xlab = "Effective sample size"
)

vioplot(
  gelman.diag(mpost_ABU$Beta, multivariate = FALSE)$psrf,
  main = "psrf(beta) - ABU model"
)

dev.off()
