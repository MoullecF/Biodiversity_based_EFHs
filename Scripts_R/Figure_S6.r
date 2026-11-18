###############################################################################
# Title: Model Fit Indicators by Life Stage (Figure S6)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())
source("Scripts_R/0-Load libraries.R")

# -----------------------------------------------------------------------------
# Load fitted models
# -----------------------------------------------------------------------------

mpa_model <- get(load("../Models/model_mpa_thin_4000_samples_250_chains_4.RData"))
abu_model <- get(load("../Models/model_abu_thin_4000_samples_250_chains_4.RData"))

# -----------------------------------------------------------------------------
# Predicted values and model fit metrics
# -----------------------------------------------------------------------------

pred_mpa <- computePredictedValues(mpa_model, expected = TRUE)
fit_mpa <- evaluateModelFit(mpa_model, predY = pred_mpa)

pred_abu <- computePredictedValues(abu_model, expected = TRUE)
fit_abu <- evaluateModelFit(abu_model, predY = pred_abu)

# -----------------------------------------------------------------------------
# Long-format performance tables
# -----------------------------------------------------------------------------

perf_mpa <- data.frame(
  species    = mpa_model$spNames,
  life_stage = substr(mpa_model$spNames, start = 9, stop = nchar(mpa_model$spNames)),
  RMSE       = fit_mpa$RMSE,
  AUC        = fit_mpa$AUC,
  TjurR2     = fit_mpa$TjurR2,
  Prevalence = colSums(mpa_model$Y) / mpa_model$ny,
  Model      = "Presence-Absence"
)

perf_mpa <- melt(perf_mpa)

perf_abu <- data.frame(
  species    = abu_model$spNames,
  life_stage = substr(abu_model$spNames, start = 9, stop = nchar(abu_model$spNames)),
  RMSE       = fit_abu$RMSE,
  R2         = fit_abu$R2,
  Prevalence = colSums(((abu_model$Y > 0) * 1) / abu_model$ny, na.rm = TRUE),
  Model      = "Conditional Abundance"
)

perf_abu <- melt(perf_abu)

performance_long <- rbind(perf_mpa, perf_abu)

# -----------------------------------------------------------------------------
# Summary tables (retained for reference/reporting)
# -----------------------------------------------------------------------------

performance_by_stage <- performance_long %>%
  dplyr::group_by(life_stage, Model, variable) %>%
  dplyr::summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n          = dplyr::n(),
    .groups    = "drop"
  )

performance_overall <- performance_long %>%
  dplyr::group_by(Model, variable) %>%
  dplyr::summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n          = dplyr::n(),
    .groups    = "drop"
  )

# -----------------------------------------------------------------------------
# Wilcoxon tests by life stage (stored for reference)
# -----------------------------------------------------------------------------

wilcox_r2 <- performance_long %>%
  filter(variable == "R2") %>%
  wilcox.test(value ~ life_stage, data = .)

wilcox_tjur <- performance_long %>%
  filter(variable == "TjurR2") %>%
  wilcox.test(value ~ life_stage, data = .)

wilcox_auc <- performance_long %>%
  filter(variable == "AUC") %>%
  wilcox.test(value ~ life_stage, data = .)

# -----------------------------------------------------------------------------
# Figure S6: Model fit indicators by life stage
# -----------------------------------------------------------------------------

metrics_to_plot <- c("AUC", "TjurR2", "R2")
performance_plot <- performance_long[performance_long$variable %in% metrics_to_plot, ]

gg_explained_power <- ggplot(performance_plot, aes(x = variable, y = value, fill = life_stage)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15),
    alpha = 0.2,
    shape = 21
  ) +
  scale_fill_manual(
    values = c("#E69F00", "#56B4E9"),
    name = "Life stage",
    labels = c("Adult", "Juvenile")
  ) +
  labs(y = "Value", x = NULL) +
  theme_bw()

# ggplot2::ggsave(gg_explained_power, filename = "../Figures/Figure_S6.png", width = 30, height = 15, units = "cm", dpi = 400)
